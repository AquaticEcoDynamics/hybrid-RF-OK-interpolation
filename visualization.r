## build the model for map2
names(M2_train)<-c("Soil", "Veg", "Landuse", "GW_depth", "Distance_GWC", "DON","s2","Latitude")

WP2Train<-reclass(as.data.frame(M2_train),a1,a2)

WP2Train<-WP2Train[,-c(4,5,10,11)]

## scale the dataset using training data 
name_list <- list("Soil", "Veg", "Landuse", "GW_depth", "Distance", "DON")
value_list <- list(b1,b2,b3)
value_max_list <-list(soil_max,veg_max,landuse_max)

for (i in seq(1,3)) {
  sub_layer <- landscapes@layers[[i]]
  sub_data <- as.data.frame(sub_layer) 
  print(i)
  for (q in seq(1, 37044)) {
    print(q)
    if ((is.numeric(sub_data[q, ])) & (sub_data[q, ] %in% value_list[[i]][, 1])) {
      sub_data[q,] <- value_list[[i]][value_list[[i]][name_list[[i]]] == sub_data[q,], 2]
    } else if ((is.numeric(sub_data[q, ])) & (!(sub_data[q, ] %in% value_list[[i]][, 1]))) {
      sub_data[q,] <- value_max_list[[i]]
    }
  }
  values(landscapes@layers[[i]]) <- sub_data[, 1]
  names(landscapes@layers[[i]]) <- name_list[[i]]
}

plot(landscapes)

WP2Train$Distance<-log10(WP2Train$Distance+0.01)

WP2Train_raster<-stack(landscapes[[1]],landscapes[[2]],landscapes[[3]],landscapes[[6]],landscapes[[7]],landscapes[[8]])
WP2Train_v<-as.data.frame(values(WP2Train_raster))
WP2Train_v$Distance<-log10(WP2Train_v$Distance+0.01)

map2_predict <- predict(rf_DON_m2, newdata = WP2Train_v)

map2<-dat.krg_DON
values(map2)<-map2_predict$data$response
#convert the raster to points for plotting
map2_df<-raster::mask(map2,study_area_withW) %>% rasterToPoints(.) %>% data.frame(.)
#Make appropriate column headings
colnames(map2_df) <- c("Longitude", "Latitude", "DON")

#Now make the map
ggplot(data = map2_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                    name = "DON mg/L",
                    breaks = c("Low", "Medium", "High"),
                    labels = c("Low", "Medium", "High"))


## map4, kriging first and then rf
# kriging for DOC
f.DOC <- as.formula(log10(DOC) ~ 1)

training_DOC <- training[,c(1,2,3,4)] %>% rbind(.,extra_n[,c(1,2,3,4)]) %>%
  rbind(.,DOC_GW4) %>% subset(.,.[,"DOC"]!="NA") %>% read_pointDataframes(.)

training_DOC<-add_S1S2(training_DOC)
var.smpl_DOC <- variogram(f.DOC, training_DOC)
plot(var.smpl_DOC)

dat.fit_DOC <- fit.variogram(var.smpl_DOC,vgm(c("Sph","Exp")))
plot(var.smpl_DOC,dat.fit_DOC)
# Perform the krige interpolation (note the use of the variogram model
dat.krg_DOC <- krige(f.DOC, training_DOC, base_grid, dat.fit_DOC) %>% raster(.) 
values(dat.krg_DOC) <- 10 ^ (values(dat.krg_DOC))

## create rasterstack with kriging data
kriging_nutrietn<-stack(dat.krg_DOC)
names(kriging_nutrietn) <- c("DOC_k")

## extract the data from landscapes_withN
landscape_train_withKN <- raster::extract(kriging_nutrietn, read_points(base6[,15:17]))

M4_train_withKN <- cbind(base6[, c(12,10,8, 6, 4,2, 13:17)],as.data.frame(landscape_train_withKN))

## create the training and testing sets 
## build the model for map2
names(M4_train_withKN)[1:11]<-c("Soil", "Veg", "Landuse","SS","GS","Catchment", "GW_depth", "Distance", "DON","Longitude","Latitude")

M4_train_withKN<-reclass(M4_train_withKN,a1,a2) %>% .[,-c(4,5,10,11)]

M4_train_withKN$DOC_SOIL<-M4_train_withKN$DOC_k*M4_train_withKN$Soil
M4_train_withKN$DOC_VEG<-M4_train_withKN$DOC_k*M4_train_withKN$Veg
M4_train_withKN$DOC_LAND<-M4_train_withKN$DOC_k*M4_train_withKN$Landuse
M4_train_withKN$DOC_CAT<-M4_train_withKN$Catchment*M4_train_withKN$DOC_k

set.seed(35)
rf_DON_m4<-model_build(M4_train_withKN,"DON","cla")

## map3 predict accuracy
map4_predict<-predict(rf_DON_m4,newdata=M4_stack_v)
#  map4_predict$data$response=map4_predict$data$response*sd_train_DON+mean_train_DON
map4<-dat.krg_DON

values(map4)<-map4_predict$data$response
#convert the raster to points for plotting
map4_df<-raster::mask(map4,study_area_withW) %>% rasterToPoints(.) %>% data.frame(.)
#Make appropriate column headings
colnames(map4_df) <- c("Longitude", "Latitude", "DON")

#Now make the map
ggplot(data = map4_df, aes(y = Latitude, x = Longitude)) +
  geom_raster(aes(fill = as.factor(DON))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#999999", "#E69F00", "#56B4E9"),
                    name = "DON mg/L",
                    breaks = c("Low", "Medium", "High"),
                    labels = c("Low", "Medium", "High")) 

### get high DON area 
high_area<-map4_df
high_area[high_area$DON==2,]<-1

ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
  geom_raster(aes(fill = as.factor(high_area$DON))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#999999","#56B4E9"),
                    name = "DON mg/L",
                    breaks = c("Low", "Medium", "High"),
                    labels = c("Low", "Medium", "High")) 

## get the uncertainty 
m4_p<-map4_predict$data
m4_p$main<-apply(m4_p[,1:3],1,max)
m4_p$uncertainty<-1-m4_p$main

map_uncertainty<-dat.krg_DON
values(map_uncertainty)<-m4_p$uncertainty
#convert the raster to points for plotting
map_un_df<-raster::mask(map_uncertainty,study_area_withW) %>% rasterToPoints(.) %>% data.frame(.)

ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
  geom_raster(aes(fill = map_un_df$var1.pred)) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  #geom_point(data=data.frame(training_df),aes(x=training_df$s1,y=training_df$s2),
  #          shape=5,size=1.5,col="red")+
  geom_point(data=data.frame(depth),aes(x=depth$s1,y=depth$s2),
             shape=3,size=2,col="red")

# map4_predict$data$response=map4_predict$data$response*(max_train_DON-min_train_DON)+min_train_DON
high_area2<-map4_df
high_area2$uncertianty <-map_un_df$var1.pred
high_area2$risk <-map_un_df$var1.pred
high_area2$needSamp <-map_un_df$var1.pred

high_area2[, "risk"][(high_area2[, "DON"] !=1)&(high_area2[,"uncertianty"]<=0.3)] <- "High_risk"
high_area2[, "risk"][high_area2[, "risk"] !="High_risk"] <- "No_risk"

high_area2[, "needSamp"][(high_area2[, "DON"] !=1)&(high_area2[,"uncertianty"]>=0.55)] <- "More_data"
high_area2[, "needSamp"][high_area2[, "needSamp"] !="More_data"] <- "No_need"

## get high DON area with low uncertatiy 
ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
  geom_raster(aes(fill = as.factor(high_area2$risk))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#56B4E9","#D0D3D4"),
                    name = "DON mg/L",
                    breaks = c("Low", "Medium", "High"),
                    labels = c("Low", "Medium", "High")) +
  geom_point(data=data.frame(depth),aes(x=depth$s1,y=depth$s2),
             shape=3,size=3,col="red")

## get median and high DON areaa with high uncertainty 

ggplot(data = map4_df, aes(y = map4_df$Latitude, x = map4_df$Longitude)) +
  geom_raster(aes(fill = as.factor(high_area2$needSamp))) + theme_bw() +
  coord_equal() +
  theme(panel.grid = element_blank(), legend.position = "right", legend.key = element_blank())+
  scale_fill_manual(values = c("#56B4E9","#D0D3D4"),
                    name = "DON mg/L",
                    breaks = c("Low", "Medium", "High"),
                    labels = c("Low", "Medium", "High"))+
  geom_point(data=data.frame(training_df),aes(x=training_df$s1,y=training_df$s2),
             shape=5,size=1)+
  geom_point(data=data.frame(depth),aes(x=depth$s1,y=depth$s2),
             shape=3,size=3,col="red")






