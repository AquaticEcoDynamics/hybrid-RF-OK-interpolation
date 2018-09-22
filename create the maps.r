
## load the library 
rm(list=ls())

library(raster)
library(sp)
library(randomForest)
library(magrittr)
library(rgdal)
library(gstat)
library(ggplot2)
library(mlr)
library(SemiPar)
library(Hmisc)
library(foreign)
library(maptools)
library(prettymapr)
library(mlrMBO)
library(parallelMap)
library(caret)
library(automap)
library(reshape2)

## start the parallel 
parallelStartSocket(4)

WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

study_area <- shapefile("study_area.shp")
water <- shapefile("water.shp")

## load the veg, soil, land use, groundwater subarea,
## surface water subare, catchment
Soil <- raster("soil1.ovr")
Veg <- raster("vegetation.ovr")
Land_use <- raster("landuse1.ovr")
Cat <- raster("catch_name2.tif.ovr")
DEM<-raster("topo_ProjectRaster3.tif.ovr")
## define the function 
## preprocess 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)

# study_area_withW
study_area_withW <- symdif(study_area, water)

pre <- function(x) {
  projection(x) <- WGS84
  extent(x) <- extent
  x <- raster::mask(x, study_area)
  return(x)
}

read_points <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SP <- spTransform(SP, WGS84)
  SP@bbox <- study_area@bbox
  if (length(zerodist(SP)) >= 1) {
    SP <- SP[-(zerodist(SP)[, 1]),]
  }
  #   plot(study_area_withW)
  #  points(SP@coords)
  return(SP)
}

read_pointDataframes <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SPD <- SpatialPointsDataFrame(SP, read_data)
  SPD <- spTransform(SPD, WGS84)
  SPD@bbox <- study_area@bbox
  if (length(zerodist(SPD)) >= 1) {
    SPD <- SPD[-(zerodist(SPD)[, 1]),]
  }
  # plot(study_area_withW)
  #points(SPD@coords)
  return(SPD)
}


reclass <- function(df, i, j) {
  df[, "DON"][df[, "DON"] <= i] <- "Low"
  df[, "DON"][df[, "DON"] < j] <- "Medium"
  df[, "DON"][(df[, "DON"] != "Low") & (df[, "DON"] != "Medium")] <- "High"
  df[, "DON"] <- factor(df[, "DON"], levels = c("Low", "Medium", "High"))
  return(df)
}


reclass4<-function(df,i,j){
  for (t in c(1,2)){
    df[, t][df[, t] <=i] <- "Low"
    df[, t][df[, t] < j] <- "Medium"
    df[, t][(df[, t] != "Low") & (df[, t] != "Medium")] <- "High"
    df[, t] <- factor(df[, t], levels = c("Low", "Medium", "High"))
  }
  return(df)
}
# Add X and Y to training 
add_S1S2 <- function(dataset) {
  dataset$s1 <- coordinates(dataset)[, 1]
  dataset$s2 <- coordinates(dataset)[, 2]
  return(dataset)
}


get_landscape<-function(df){
  landscape_all<-data.frame()
  for (ii in seq(1,length(df))){
    aa<-as.data.frame(df[[ii]])
    aa<-subset(aa,aa$Soil!="NA")
    Soil=tail(names(sort(table(aa[,1]))),1)
    Veg=tail(names(sort(table(aa[,2]))),1)
    Landuse=tail(names(sort(table(aa[,3]))),1)
    Catchment=tail(names(sort(table(aa[,4]))),1)
    GW_depth=mean(aa[,5])
    Distance=mean(aa[,6])
    Distance_LP=mean(aa[,7])
    Distance_GWC=mean(aa[,8])
    slope=mean(aa[,9])
    aspect=mean(aa[,10])
    sing_land<-data.frame(Soil,Veg,Landuse,Catchment,GW_depth,Distance,Distance_LP,Distance_GWC,slope,aspect)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}



get_landscape2<-function(df){
  landscape_a<-data.frame()
  for (ii in seq(1,length(df))){
    aa<-as.data.frame(df[[ii]])
    aa<-subset(aa,aa$Soil!="NA")
    Soil=tail(names(sort(table(aa[,1]))),1)
    Veg=tail(names(sort(table(aa[,2]))),1)
    Landuse=tail(names(sort(table(aa[,3]))),1)
    GW_depth=mean(aa[,4])
    Distance_GWC=mean(aa[,5])
    slope=mean(aa[,6])
    s1=mean(aa[,7])
    s2=mean(aa[,8])
    sing_land<-data.frame(Soil,Veg,Landuse,GW_depth,Distance_GWC,slope,s1,s2)
    landscape_a<-rbind(landscape_a,sing_land)
  }
  return(landscape_a)
}

## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)
DEM <- pre(DEM)

v_Veg<-values(Veg)
v_Veg[v_Veg %in% c(2,3,4)]=1
v_Veg[v_Veg %in% c(8,9)]=8
v_Veg[v_Veg %in% c(12,13)]=12
v_Veg[v_Veg %in% c(18,19,20)]=18
values(Veg)<-v_Veg

v_land<-values(Land_use)
v_land[v_land %in% c(1,2,13)]=1
v_land[v_land %in% c(5,7,12,6,11)]=5
v_land[v_land %in% c(8,10)]=8
values(Land_use)<-v_land

v_soil<-values(Soil)
v_soil[v_soil %in% c(1,2)]=1
v_soil[v_soil %in% c(4,5)]=4
v_soil[v_soil %in% c(6,7)]=6
v_soil[v_soil %in% c(11,12)]=11
v_soil[v_soil %in% c(13,14)]=13
values(Soil)<-v_soil

# Create an empty grid where n is the total number of cells
r <- raster(study_area)
res(r) <- res(Soil) # 10 km if your CRS's units are in km
base_grid <- as(r, 'SpatialGrid')
#plot(base_grid)

## M2, using RF to predict the DON
depth <- read.csv("sampling_depth.csv",header=T) %>% read_pointDataframes(.)

# Define the 1st order polynomial equation
f_depth <- as.formula(sampling_d ~ 1)
# Add X and Y to training 
depth<-add_S1S2(depth)
# variogram on the de-trended data.
var.depth <- variogram(f_depth, depth)
plot(var.depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp")))
plot(var.depth,dat.fit_depth)

# created in the earlier step)
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.)
#plot(depth_k)
depth_k@data@names<-"GW_depth"

#Now make the map
### distance
water <- raster::rasterize(water, depth_k)
#water_distance <- raster::mask(distance(water),study_area)
water_distance <- raster::distance(water)

water_distance@data@names<-"Distance_to_water"

GW_center<-data.frame(Latitude=c(6495000,6475000,6460000,6448000,6403000),Longitude=rep(402000,5),values=1)
GW_center <- SpatialPoints(GW_center[, c(2:1)], proj4string = WGS84)
GW_center@bbox <- study_area@bbox
base_GWC<-water 
values(base_GWC)<-1
Distance_GWC<-distanceFromPoints(base_GWC,GW_center)
Distance_GWC@data@names<-"Distance_GWC"

## load the data 
landscapes<-stack(Soil,Veg,Land_use,depth_k,water_distance,Distance_GWC)
names(landscapes) <- c("Soil", "Veg", "Landuse", "GW_depth", "Distance","Distance_GWC")

## load the data 
set.seed(666)

all_points<-read.csv("all_data1127.csv",header = T)
extra_n<-read.csv("extra_n.csv",header = T)
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))
extra_n<-subset(extra_n,!(extra_n$WIN_Site_ID %in% all_points$WIN_Site_ID))

a1=0.5
a2=2.0

training<-all_points

#training<-all_data[ trainIndex,]

## load the point data 
training_df <- read_pointDataframes(training)
#testing_df <-  read_pointDataframes(testing) 

training_points<- read_points(training)
#testing_points <- read_points(testing)

## map1, using kringing for DON interpolation
f.1 <- as.formula(log(DON) ~ 1)
# Add X and Y to training 
training_df<-add_S1S2(training_df)
#testing_df<-add_S1S2(testing_df)

# Compute the sample variogram; note that the f.1 trend model is one of the
var.smpl1 <- variogram(f.1, training_df)
# plot(var.smpl1)
# Compute the variogram model by passing the nugget, sill and range value
dat.fit1 <-   fit.variogram(var.smpl1,fit.sills = FALSE, fit.ranges = FALSE, 
                            model = vgm(nugget = 0.5, "Exp", range = 3000,  psill = 0.45)) 

#options(repr.plot.width=5, repr.plot.height=5)  
plot(var.smpl1,dat.fit1)

# Perform the krige interpolation (note the use of the variogram model
kriging_DON_m1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) 
values(kriging_DON_m1)<-10^(values(kriging_DON_m1))

dat.krg_DON<-kriging_DON_m1

#convert the raster to points for plotting
DON_map1_2<-raster::mask(kriging_DON_m1,study_area_withW)  
map1_df <- rasterToPoints(DON_map1_2) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(map1_df) <- c("Longitude", "Latitude", "DON")
t=3
map1_df[, t][map1_df[, t] <=a1] <- "Low"
map1_df[, t][map1_df[, t] < a2] <- "Medium"
map1_df[, t][(map1_df[, t] != "Low") & (map1_df[, t] != "Medium")] <- "High"
map1_df[, t] <- factor(map1_df[, t], levels = c("Low", "Medium", "High"))


#Now make the map
p11<-ggplot(data = map1_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#B8B8B8", "#FFD700", "#FF3030"),
                    name = "DON concentration",
                    breaks = c( "High","Medium", "Low"),
                    labels = c("High","Medium", "Low"))+labs(x=as.character(a1))
plot(p11)

#convert the raster to points for plotting
map1_df2 <- rasterToPoints(DON_map1_2) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(map1_df2) <- c("Longitude", "Latitude", "DON")

#Now make the map
ggplot(data = map1_df2, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = DON)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
scale_fill_gradientn(colours = c("#B8B8B8", "#FFD700", "#FF3030"),breaks=c(0,1,2,3,4,5,6,7),limits=c(0,7))


landscape_train <- capture_zone_land(training_df)

M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON")])

common_landscape<-function(land){
  land_dataset<-data.frame(table(M2_train[,land]))
  land_common<-subset(land_dataset,land_dataset[,2]==max(land_dataset[,2]))[1]
  return(as.matrix(land_common))
}

soil_max = common_landscape("Soil")[1]
veg_max=common_landscape("Veg")[1]
landuse_max = common_landscape("Landuse")[1]

max_list<-list(soil_max,veg_max,landuse_max)

b1=data.frame(table(M2_train$Soil))
b2=data.frame(table(M2_train$Veg))
b3=data.frame(table(M2_train$Landuse))
value_list<-list(b1,b2,b3)
name_list <- list("Soil", "Veg", "Landuse")

landscapes_all<-landscapes
for (ii in seq(1,3,1)) {
  sub_layer <- landscapes_all@layers[[ii]]
  sub_data <- as.data.frame(sub_layer) 
  print(ii)
  for (q in seq(1, 37044)) {
    if ((is.numeric(sub_data[q, ])) & (!(sub_data[q, ] %in% value_list[[ii]][, 1]))){
      sub_data[q,1] <- as.numeric(max_list[[ii]])
      print(q)
    }
  }
  values(landscapes_all@layers[[ii]]) <- sub_data[, 1]
  names(landscapes_all@layers[[ii]]) <- name_list[[ii]]
}

plot(landscapes_all)

M2_map<-stack(landscapes_all[[1]],landscapes_all[[2]],landscapes_all[[3]],landscapes_all[[4]],
              landscapes_all[[5]],landscapes_all[[6]])

#WP2Train_v<-as.data.frame(values(M2_map))

#map2<-raster::mask(landscapes[[1]],study_area_withW)
#convert the raster to points for plotting
map2_all_points<- rasterToPoints(M2_map) %>% data.frame(.) 
map2_all_points<-map2_all_points[,c(1,2,3)]
#Make appropriate column headings
colnames(map2_all_points) <- c("Longitude", "Latitude", "DON")
SP <- SpatialPoints(map2_all_points[, 1:2], proj4string =  WGS84)
SPD <- SpatialPointsDataFrame(SP, map2_all_points)
SPD@bbox <- study_area@bbox


a=100
b=200
capture_zone_land2<-function(df){
  num<-nrow(df)
  landscape_data<-data.frame()
  for (r in seq(1,num)){
    print(r)
    p1_long<-df@coords[r,1]
    p1_lat<-df@coords[r,2]
    pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long+a,p1_lat+b),c(p1_long+2*a,p1_lat+b),
                         c(p1_long+2*a,p1_lat-b),c(p1_long+a,p1_lat-b),c(p1_long,p1_lat)))  
    projection(pg)<- WGS84
    p1_landscape<-raster::extract(M2_map,pg)
    p1_landscape<-get_landscape(p1_landscape)
    landscape_data<-rbind(landscape_data,p1_landscape)
  }
  return(landscape_data)
}


landscape_all_points <- capture_zone_land2(SPD)

landscape_all_points2<-cbind(landscape_all_points)

M2_map_data<-landscape_all_points2

WP2Train<-reclass(M2_train,a1,a2)

WP2Train<- as.h2o(WP2Train)
M2_map_data<- as.h2o(M2_map_data)

y <- "DON"
x <- setdiff(names(WP2Train), y)

rf_m2 <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = WP2Train,
                          seed = 123)  

map2_predict=as.data.frame(h2o::h2o.predict(rf_m2,M2_map_data))

map2<-landscapes_all[[1]]
values(map2)<-map2_predict$predict

DON_map2_2<-raster::mask(map2,study_area_withW)

map2_df <- rasterToPoints(DON_map2_2) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(map2_df) <- c("Longitude", "Latitude", "DON")
t=3
map2_df[, t][map2_df[, t] ==1] <- "High"
map2_df[, t][map2_df[, t] ==2 ] <- "Low"
map2_df[, t][(map2_df[, t] != "High") & (map2_df[, t] != "Low")] <- "Medium"
map2_df[, t] <- factor(map2_df[, t], levels = c("Low", "Medium", "High"))

#Now make the map
ggplot(data = map2_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#B8B8B8", "#FFD700", "#FF3030"),
                    name = "DON concentration",
                    breaks = c( "High","Medium", "Low"),
                    labels = c("High","Medium", "Low"))


## kriging
f_DOC <- as.formula(log10(DOC) ~ 1)  
var_DOC <- variogram(f_DOC, training_df)

dat_DOC <-   fit.variogram(var_DOC,fit.sills = FALSE, fit.ranges = FALSE, 
                           model = vgm(nugget = 0.15, "Exp", range = 15000,  psill = 0.4))

#options(repr.plot.width=5, repr.plot.height=5)  
plot(var_DOC,dat_DOC)

kriging_DOC <- krige(f_DOC, training_df, base_grid, dat_DOC) %>% raster(.) 
values(kriging_DOC) <- 10 ^ (values(kriging_DOC))

## create rasterstack with kriging data
kriging_nutrietn_DOC<-stack(kriging_DOC,kriging_DON_m1)
names(kriging_nutrietn_DOC) <- c("DOC_k","DON_k")


## extract the data from landscapes_withN
landscape_train_withKN <- raster::extract(kriging_nutrietn_DOC,training_df)
landscape_all_points_DOC <-  raster::extract(kriging_nutrietn_DOC,SPD)


M4_train_withKN<- cbind(as.data.frame(WP2Train),as.data.frame(landscape_train_withKN)) 
M4_train_all <- cbind(as.data.frame(M2_map_data),as.data.frame(landscape_all_points_DOC)) 

y <- "DON"
x <- setdiff(names(M4_train_withKN), y)

M4_train_withKN<-as.h2o(M4_train_withKN)
M4_train_all<-as.h2o(M4_train_all)

# GBM hyperparamters
rf_m4 <- h2o.randomForest(x = x,
                          y = y,
                          #nfolds = 3,
                          training_frame = M4_train_withKN,
                          seed = 123)  

# Now let's evaluate the model performance on a test set

map4_predict=as.data.frame(h2o::h2o.predict(rf_m4,M4_train_all))

map4<-landscapes_all[[1]]
values(map4)<-map4_predict$predict

DON_map4<-raster::mask(map4,study_area_withW)
map4_df <- rasterToPoints(DON_map4) %>% data.frame(.)

#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(map4_df) <- c("Longitude", "Latitude", "DON")
t=3
map4_df[, t][map4_df[, t] ==1] <- "High"
map4_df[, t][map4_df[, t] ==2 ] <- "Low"
map4_df[, t][(map4_df[, t] != "High") & (map4_df[, t] != "Low")] <- "Medium"
map4_df[, t] <- factor(map4_df[, t], levels = c("Low", "Medium", "High"))


#Make the points a dataframe for ggplot
#Now make the map
ggplot(data = map4_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#B8B8B8", "#FFD700", "#FF3030"),
                    name = "DON concentration",
                    breaks = c( "High","Medium", "Low"),
                    labels = c("High","Medium", "Low"))



GW_depth_df <- raster::mask(landscapes[[4]],study_area_withW) %>% 
  rasterToPoints(.) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(GW_depth_df) <- c("Longitude", "Latitude", "GW_depth")
GW_depth_df$GW_depth<--GW_depth_df$GW_depth
#Now make the map
ggplot(data = GW_depth_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = GW_depth)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_continuous(low = "red", high = "yellow")


DOC_kriging<- raster::mask(kriging_DOC,study_area_withW) %>% 
  rasterToPoints(.) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(DOC_kriging) <- c("Longitude", "Latitude", "DOC")

#Now make the map
ggplot(data = DOC_kriging, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = DOC)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_continuous(low = "red", high = "yellow")


Distance_df <- raster::mask(landscapes[[5]],study_area_withW) %>% 
  rasterToPoints(.) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(Distance_df) <- c("Longitude", "Latitude", "distance")

#Now make the map
ggplot(data = Distance_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = distance)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_gradientn(colours = terrain.colors(5),limits=c(0,60000))


Distance_gwc_df <- raster::mask(landscapes[[6]],study_area_withW) %>% 
  rasterToPoints(.) %>% data.frame(.)
#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(Distance_gwc_df) <- c("Longitude", "Latitude", "distance_gwc")

#Now make the map
ggplot(data = Distance_gwc_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = distance_gwc)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_gradientn(colours = terrain.colors(5),limits=c(0,60000))


### uncertainty map

#map4_predict<-predict(rf_DON_m4,newdata=M4_train_all)
map4_uncertainty<-landscapes_all[[1]]

map4_predict$main<- apply(map4_predict[, 2:4], 1, max)
map4_predict$uncertainty<- 1-map4_predict$main

values(map4_uncertainty)<-map4_predict$uncertainty

uncer_map4<-raster::mask(map4_uncertainty,study_area_withW)
map4_uncert_df <- rasterToPoints(uncer_map4) %>% data.frame(.)

#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(map4_uncert_df) <- c("Longitude", "Latitude", "Uncertainty")

#Make the points a dataframe for ggplot
#Now make the map
training_dataset<-as.data.frame(training_df@data)

ggplot(data = map4_uncert_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = Uncertainty)) + theme_bw() +
  geom_point(data=training_dataset,aes(y=training_dataset$s2,x=training_dataset$s1),
             col="black",shape=3,size=2.5)+
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right")+
  scale_fill_continuous(low = "#AED6F1", high = "#21618C")

## uncertainty map_2
map2_uncertainty<-landscapes_all[[1]]

map2_predict$main<- apply(map2_predict[, 2:4], 1, max)
map2_predict$uncertainty<- 1-map2_predict$main

values(map2_uncertainty)<-map2_predict$uncertainty

uncer_map2<-raster::mask(map2_uncertainty,study_area_withW)
map2_uncert_df <- rasterToPoints(uncer_map2) %>% data.frame(.)

#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(map2_uncert_df) <- c("Longitude", "Latitude", "Uncertainty")

#Make the points a dataframe for ggplot
#Now make the map
ggplot(data = map2_uncert_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = Uncertainty)) + theme_bw() +
  geom_point(data=training_dataset,aes(y=training_dataset$s2,x=training_dataset$s1),
             size=2.5,shape=3,col="black")+
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_continuous(low = "#AED6F1", high = "#21618C")


### high risk area
high_risk<-landscapes_all[[1]]

map4_predict$high_risk <- ifelse(((map4_predict$predict=="High")&(map4_predict$uncertainty)<0.3),
                                 c("High risk"), c("Other")) 

map4_predict$high_risk<-factor(map4_predict$high_risk)

values(high_risk)<-map4_predict$high_risk

high_risk<-raster::mask(high_risk,study_area_withW)
high_risk <- rasterToPoints(high_risk) %>% data.frame(.)

#Make the points a dataframe for ggplot
#Make appropriate column headings
colnames(high_risk) <- c("Longitude", "Latitude", "Risk")
t=3
high_risk[, t][high_risk[, t] ==1] <- "High risk"
high_risk[, t][(high_risk[, t] != "High risk") ] <- "Other areas"
high_risk[, t] <- factor(high_risk[, t], levels = c("High risk", "Other areas"))


#Make the points a dataframe for ggplot
#Now make the map

ggplot(data = high_risk, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = as.factor(Risk))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#FF3030","#B8B8B8"),
                    name = "DON concentration",
                    breaks = c( "High risk", "Other areas"),
                    labels =  c( "High risk", "Other areas"))

