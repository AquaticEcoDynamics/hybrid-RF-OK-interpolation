
## load the library 
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
library(rminer)

## start the parallel 
h2o.init(nthreads = 16,max_mem_size = "40g")

WGS84 <- CRS("+proj=utm +zone=50 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
study_area <- shapefile("study_area.shp")
water <- shapefile("water.shp")

## load the veg, soil, land use
Soil <- raster("soil1.ovr")
Veg <- raster("vegetation.ovr")
Land_use <- raster("landuse1.ovr")
Cat <- raster("catch_name2.tif.ovr")
DEM<-raster("topo_ProjectRaster3.tif.ovr")


## define the function 
study_area <- spTransform(study_area, WGS84)
extent <- c(study_area@bbox[1, 1:2], study_area@bbox[2, 1:2])

water <- spTransform(water, WGS84)
## preprocess 
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
    GW_depth=mean(aa[,4])
    Distance=mean(aa[,5])
    Distance_GWC=mean(aa[,6])
    sing_land<-data.frame(Soil,Veg,Landuse,GW_depth,Distance,Distance_GWC)
    landscape_all<-rbind(landscape_all,sing_land)
  }
  return(landscape_all)
}

read_points <- function(read_data) {
  SP <- SpatialPoints(read_data[, 2:3], proj4string = CRS("+proj=longlat +ellps=GRS80 +no_defs"))
  SP <- spTransform(SP, WGS84)
  SP@bbox <- study_area@bbox
  if (length(zerodist(SP)) >= 1) {
    SP <- SP[-(zerodist(SP)[, 1]),]
  }
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
  return(SPD)
}

## M2, using RF to predict the DON
a=100
b=200
capture_zone_land<-function(df){
  num<-nrow(df)
  landscape_data<-data.frame()
  for (r in seq(1,num)){
    print(r)
    p1_long<-df@coords[r,1]
    p1_lat<-df@coords[r,2]
    pg<-spPolygons(rbind(c(p1_long,p1_lat),c(p1_long+a,p1_lat+b),c(p1_long+2*a,p1_lat+b),
                         c(p1_long+2*a,p1_lat-b),c(p1_long+a,p1_lat-b),c(p1_long,p1_lat)))  
    projection(pg)<- WGS84
    p1_landscape<-raster::extract(landscapes,pg)
    p1_landscape<-get_landscape(p1_landscape)
    landscape_data<-rbind(landscape_data,p1_landscape)
  }
  return(landscape_data)
}

## preprocess the landscape raster
Soil <- pre(Soil)
Veg <- pre(Veg)
Land_use <- pre(Land_use)

v_Veg<-values(Veg)
v_Veg[v_Veg %in% c(2,3,4)]=1
v_Veg[v_Veg %in% c(8,9)]=8
v_Veg[v_Veg %in% c(12,13)]=12
v_Veg[v_Veg %in% c(18,19,20)]=18
values(Veg)<-v_Veg

v_land<-values(Land_use)
v_land[v_land %in% c(1,2,5,6,7,11,12,13)]=1
v_land[v_land %in% c(3,4)]=3
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
res(r) <- res(Soil) 
base_grid <- as(r, 'SpatialGrid')

## M2, using RF to predict the DON
depth <- read.csv("depth_to_groundwater.csv",header=T) %>% read_pointDataframes(.)

# Define the 1st order polynomial equation
f_depth <- as.formula(sampling_d ~ 1)
# Add X and Y to training 
depth<-add_S1S2(depth)
# variogram on the de-trended data.
var.depth <- variogram(f_depth, depth)
#plot(var.depth)
dat.fit_depth <- fit.variogram(var.depth,vgm(c("Sph","Exp")))
# created in the earlier step)
depth_k <- krige(f_depth, depth, base_grid, dat.fit_depth) %>% raster(.) %>% raster::mask(., study_area)
#plot(depth_k)
depth_k@data@names<-"GW_depth"

#Now make the map
### distance
water <- raster::rasterize(water, depth_k)
water_distance <- raster::mask(distance(water),study_area)
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

## 
## load the data 
set.seed(666)

all_results<-data.frame()
all_data<-read.csv("groundwater_nutrient.csv",header = T)

results<-data.frame()
set.seed(91)
seed.list<-sample(1:1000,300,replace =F)

a1=0.5
a2=2.0

set.seed(91)
trainIndex <- createDataPartition(all_data$DON, p = .75, list = FALSE,time=1)

training<-all_data[ trainIndex,]
testing<-all_data[-trainIndex,]

## load the point data 
training_df <- read_pointDataframes(training)
testing_df <-  read_pointDataframes(testing) 

training_points<- read_points(training)
testing_points <- read_points(testing)

## map1, using kringing for DON interpolation
f.1 <- as.formula(log(DON) ~ 1)
# Add X and Y to training 
training_df<-add_S1S2(training_df)
testing_df<-add_S1S2(testing_df)

# Compute the sample variogram; note that the f.1 trend model is one of the
var.smpl1 <- variogram(f.1, training_df)
# Compute the variogram model by passing the nugget, sill and range value
dat.fit1 <-   fit.variogram(var.smpl1,fit.sills = FALSE, fit.ranges = FALSE, 
                            model = vgm(nugget = 0.55, "Exp", range = 15000,  psill = 0.4)) 

# Perform the krige interpolation (note the use of the variogram model
kriging_DON_m1 <- krige(f.1, training_df, base_grid, dat.fit1) %>% raster(.) %>% raster::mask(., study_area)
values(kriging_DON_m1)<-10^(values(kriging_DON_m1))

dat.krg_DON<-kriging_DON_m1

map1_predict_train <- data.frame(observed_DON=training_df@data$DON,predicted_DON=raster::extract(kriging_DON_m1, training_points))
map1_predict_train<-reclass4(map1_predict_train,a1,a2)

map1_predict <- data.frame(observed_DON=testing_df@data$DON,predicted_DON=raster::extract(kriging_DON_m1, testing_points))  
map1_predict<-reclass4(map1_predict,a1,a2)

M1_ACC<-postResample(map1_predict[,2],map1_predict[,1])[1]
M1_kappa<-postResample(map1_predict[,2],map1_predict[,1])[2]  

landscape_train <- capture_zone_land(training_df)
landscape_test <- capture_zone_land(testing_df)

M2_train <- cbind(as.data.frame(landscape_train), training_df@data[c("DON")])
M2_test <- cbind(as.data.frame(landscape_test), testing_df@data[c("DON")])

names(M2_train) <- colnames(M2_test)

common_landscape<-function(land){
  land_dataset<-data.frame(table(M2_train[,land]))
  land_common<-subset(land_dataset,land_dataset[,2]==max(land_dataset[,2]))[1]
  return(as.matrix(land_common))
}

soil_max = common_landscape("Soil")[1]
veg_max=common_landscape("Veg")[1]
landuse_max = common_landscape("Landuse")[1]

max_list<-list(soil_max,veg_max,landuse_max)

for (ii in seq(1,3)){
  M2_train[,ii]<-factor(M2_train[,ii],levels = unique(values(landscapes[[ii]]))[-1])
  M2_test[,ii]<-factor(M2_test[,ii],levels=unique(values(landscapes[[ii]]))[-1])
  M2_test [(which(!(M2_test[,ii] %in% M2_train[,ii]))),ii]<-as.numeric(max_list[[ii]])
  
  M2_train[,ii]<-droplevels(M2_train[,ii])
  M2_test[,ii]<-factor(M2_test[,ii],levels = levels(M2_train[,ii]))
}

M2_train<-reclass(M2_train,a1,a2)
M2_test<-reclass(M2_test,a1,a2)

WP2Train<- as.h2o(M2_train)
WP2Test<- as.h2o(M2_test)

y <- "DON"
x <- setdiff(names(WP2Train), y)

rf_m2 <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = WP2Train,          
                          seed = 123)  

h20_rf_pred=as.data.frame(h2o::h2o.predict(rf_m2,WP2Test))
h20_rf_pred_train=as.data.frame(h2o::h2o.predict(rf_m2,WP2Train))

WP2Test=as.data.frame(WP2Test)

M2_ACC<-postResample(h20_rf_pred$predict,WP2Test$DON)[1]
M2_kappa<-postResample(h20_rf_pred$predict,WP2Test$DON)[2]

## map1, using kringing for DOC interpolation
f_DOC <- as.formula(log10(DOC) ~ 1)  
var_DOC <- variogram(f_DOC, training_df)

dat_DOC <-   fit.variogram(var_DOC,fit.sills = FALSE, fit.ranges = FALSE, 
                           model = vgm(nugget = 0.2, "Exp", range = 15000,  psill = 0.4))

kriging_DOC <- krige(f_DOC, training_df, base_grid, dat_DOC) %>% raster(.) %>% raster::mask(., study_area)
values(kriging_DOC) <- 10 ^ (values(kriging_DOC))

kriging_DON<-stack(kriging_DOC,kriging_DON_m1)
names(kriging_DON) <- c("DOC_k","DON_k")
## extract the data from landscapes
landscape_train_withKN <- raster::extract(kriging_DON,training_df)
landscape_test_withKN <-  raster::extract(kriging_DON,testing_df)

M4_train<- cbind(as.data.frame(M2_train),as.data.frame(landscape_train_withKN)) 
M4_test <- cbind(as.data.frame(M2_test),as.data.frame(landscape_test_withKN)) 

names(M4_train)<-names(M4_test)

y <- "DON"
x <- setdiff(names(M4_train), y)

M4_train<-as.h2o(M4_train)
M4_test<-as.h2o(M4_test)

rf_m4 <- h2o.randomForest(x = x,
                          y = y,
                          training_frame = M4_train,
                          seed = 123)  

# Now let's evaluate the model performance on a test set
h20_rf_pred_m4=as.data.frame(h2o::h2o.predict(rf_m4,M4_test))

M4_test=as.data.frame(M4_test)

M4_ACC<-postResample(h20_rf_pred_m4$predict,M4_test$DON)[1]
M4_kappa<-postResample(h20_rf_pred_m4$predict,M4_test$DON)[2]

sing_acc<-data.frame(M1_ACC,M2_ACC,M4_ACC,M1_kappa,M2_kappa,M4_kappa)
results<-rbind(results,sing_acc)
print(results)
