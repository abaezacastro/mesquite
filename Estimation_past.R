require("raster")
metodo="GLM"
load("data/raster_Tmean_historic")
load("data/annual_total_rainfall_raster_historico")
load("data/raster_Tmin_historic")
multiple_rainfall=PPT_m_H*PPT_m_H
##########################################################################
#re-clasification
m <- c(-1, 0, NA,  0.1, 5000, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(PPT_m_H, rclmat)

total_annual_rainfall=PPT_m_H*rc
raster_Tmean_annual_H=raster_Tmean_annual_H*rc
multiple_rainfall=multiple_rainfall*rc
raster_Tmin_annual_H=raster_Tmin_annual_H*rc
##########################################################################
load("data/Prosopis_glandulosa_presence_data")
##########################################################################

##########################################################################
predictors_past<-stack(raster_Tmean_annual_H,
                  raster_Tmin_annual_H,
                  total_annual_rainfall,
                  multiple_rainfall
                  )

predictors_past_wo_tmin<-stack(raster_Tmean_annual_H,
                          total_annual_rainfall,
                  multiple_rainfall
)

rm(raster_Tmean_annual)
rm(raster_Tmin_annual)
rm(PPT_m_H)
rm(multiple_rainfall)
rm(total_annual_rainfall)

llmespg<-cbind(mesquite_pg_clean$lon[which(mesquite_pg_clean$year<=1950)], mesquite_pg_clean$lat[which(mesquite_pg_clean$year<=1950)])
llmespg<-llmespg[,-which(llmespg[,2]<0)]
presvals_pg<-as.data.frame(extract(predictors_past, llmespg))
presvals_pg_wo_tmin<-as.data.frame(extract(predictors_past_wo_tmin, llmespg))
names(presvals_pg)<-c("mean_annual_temperature","min_annual_temperature","total_annual_rainfall","multiple_rainfall")
names(predictors_past)<-names(presvals_pg)

names(presvals_pg_wo_tmin)<-c("mean_annual_temperature","total_annual_rainfall","multiple_rainfall")
names(predictors_past_wo_tmin)<-names(presvals_pg_wo_tmin)


##########################################################################

#training set (absences)
#backgr <- dismo::randomPoints(predictors_past, 2000)
backgr<-readRDS("data/background_points")
absvals <- as.data.frame(extract(predictors_past, backgr))
names(absvals)<-names(presvals_pg)

pg <- c(rep(1, nrow(presvals_pg)), rep(0, nrow(absvals)))
sdmdata_pg <- data.frame(cbind(pg, rbind(presvals_pg, absvals)))

##########################################################################
valid_points=which(sdmdata_pg$mean_annual_temperature!="NaN")
traning<-sample(x = (1:length(sdmdata_pg$mean_annual_temperature))[valid_points],size =round(length(valid_points)/2),replace = F)

traning_point=sdmdata_pg[traning,]
testing_points=sdmdata_pg[-traning,]

p_training=sdmdata_pg$pg[traning]



require(MASS)
#Estimations and predictions
if (metodo=="GLM"){
m1 <- glm(pg ~ mean_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m1_predict=predict(object = m1,newdata=as.data.frame(predictors_past))

m2 <- glm(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m2_predict=predict(object = m2,newdata=as.data.frame(predictors_past))
}

if (metodo=="GAM"){
m3<-gam::gam(pg ~ mean_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m3_predict<-gam::predict.Gam(object = m3,newdata=as.data.frame(predictors_past))

m4<-gam::gam(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m4_predict<-gam::predict.Gam(object = m4,newdata=as.data.frame(predictors_past))

}


if (metodo=="BIOCLIM"){
###########################################################
#bioclim uses only- presence data
valid_points_wo_tmin=which(presvals_pg_wo_tmin$mean_annual_temperature!="NaN")
traning_wo_tmin<-sample(x = (1:length(presvals_pg_wo_tmin$mean_annual_temperature))[valid_points_wo_tmin],size =round(length(valid_points_wo_tmin)/2),replace = F)
traning_point_wo_tmin=presvals_pg_wo_tmin[traning_wo_tmin,]
testing_points_wo_tmin=presvals_pg_wo_tmin[-traning_wo_tmin,]

m5<-dismo::bioclim(traning_point_wo_tmin)
m5_predict<-predict(m5,predictors_past_wo_tmin)



valid_points=which(presvals_pg$mean_annual_temperature!="NaN")
traning<-sample(x = (1:length(presvals_pg$mean_annual_temperature))[valid_points],size =round(length(valid_points)/2),replace = F)
traning_point=presvals_pg[traning,]
testing_points=presvals_pg[-traning,]

m6<-dismo::bioclim(traning_point)
m6_predict<-predict(m6,predictors_past)
}
##########################################################################



#evaluation
nn<-raster(x = predictors_past)
if (metodo=="GLM"){
dist_m1<-setValues(x = nn,values = m1_predict)
dist_m2<-setValues(x = nn,values = m2_predict)
}
if (metodo=="GAM"){
dist_m3<-setValues(x = nn,values = m3_predict)
dist_m4<-setValues(x = nn,values = m4_predict)
}
if (metodo=="BIOCLIM"){
  dist_m5<-m5_predict
dist_m6<-m6_predict
}

##############################################################################################################
USA_map_full<-maptools::readShapeSpatial("C:/Users/abaezaca/Dropbox (ASU)/Design Principles in CIS/maps_USA/tl_2017_us_state/tl_2017_us_state")
top = 49.3457868 # north lat
left = -124.7844079 # west long
right = -66.9513812 # east long
bottom =  24.7433195 # south lat
USA_map<-crop(USA_map_full,raster::extent(left,right,bottom,top))
USA_map@data$NAME<-factor(USA_map@data$NAME)
##############################################################################################################

if (metodo=="GLM"){
  dist_m1=mask(dist_m1,USA_map)
  dist_m2=mask(dist_m2,USA_map)
}
if (metodo=="GAM"){
  dist_m3=mask(dist_m3,USA_map)
  dist_m4=mask(dist_m4,USA_map)
}
if (metodo=="BIOCLIM"){
  dist_m5=mask(dist_m5,USA_map)
dist_m6=mask(dist_m6,USA_map)
}

if (metodo=="GLM"){
eval_m1<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m1,x=predictors_past)
eval_m2<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m2,x=predictors_past)
}
if (metodo=="GAM"){
eval_m3<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m3,x=predictors_past)
eval_m4<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m4,x=predictors_past)
}

if (metodo=="BIOCLIM"){
  eval_m5<-dismo::evaluate(p=llmespg[-traning_wo_tmin,],a=backgr,model = m5,x=predictors_past_wo_tmin)
eval_m6<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m6,x=predictors_past)
}


source("future_scenarios.R")
nn_future<-raster(x = predictors_future)

if (metodo=="GLM"){
  T50<-dist_m2>eval_m2@t[which.max(eval_m2@TPR+eval_m2@TNR)]
  m2_predict_T2000=predict(object = m2,newdata=as.data.frame(predictors_current))
  dist_m2<-setValues(x = nn,values = m2_predict_T2000)
  T2000<-dist_m2>eval_m2@t[which.max(eval_m2@TPR+eval_m2@TNR)]
  
  m2_predict_future_T2100=predict(object = m2,newdata=as.data.frame(predictors_future))
  dist_m2_future<-setValues(x = nn_future,values = m2_predict_future_T2100)
  dist_m2_future=mask(dist_m2_future,USA_map)
  T2100<-dist_m2_future>eval_m2@t[which.max(eval_m2@TPR+eval_m2@TNR)]
}


###############################################################################3
# png(filename = "output_maps/comparison_withTmin_past.png",width = 20,height = 30,res=300,units = "cm")
# par(mfrow=c(3,2))
# plot(dist_m1>eval_m1@t[which.max(eval_m1@TPR+eval_m1@TNR)],col=c("grey","red"),main="GLM without Tmin")
# plot(USA_map,add=T)
# plot(dist_m2>eval_m2@t[which.max(eval_m2@TPR+eval_m2@TNR)],col=c("grey","red"),main="GLM with Tmin")
# plot(USA_map,add=T)
# plot(dist_m3>eval_m3@t[which.max(eval_m3@TPR+eval_m3@TNR)],col=c("grey","red"),main="GAM without Tmin")
# plot(USA_map,add=T)
# plot(dist_m4>eval_m4@t[which.max(eval_m4@TPR+eval_m4@TNR)],col=c("grey","red"),main="GAM with Tmin")
# plot(USA_map,add=T)
# plot(dist_m5>eval_m5@t[which.max(eval_m5@TPR+eval_m5@TNR)],col=c("grey","red"),main="GAM without Tmin")
# plot(USA_map,add=T)
# plot(dist_m6>eval_m6@t[which.max(eval_m6@TPR+eval_m6@TNR)],col=c("grey","red"),main="GAM with Tmin")
# plot(USA_map,add=T)
# 
# dev.off()


##############################################################################

