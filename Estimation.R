load("data/raster_Tmean_current")
load("data/annual_total_rainfall_raster_current")
load("data/raster_Tmin_current")
multiple_rainfall=raster_annualPPT*raster_annualPPT
raster_Tmin_annual<-calc(raster_Tmin_annual,fun=function(x){x/10}) #K to C
raster_Tmean_annual<-calc(raster_Tmean_annual,fun=function(x){x/10}) #K to C
##########################################################################
#re-clasification
m <- c(-1, 0, NA,  0.1, 5000, 1)
rclmat <- matrix(m, ncol=3, byrow=TRUE)
rc <- reclassify(raster_annualPPT, rclmat)

total_annual_rainfall=raster_annualPPT*rc
raster_Tmean_annual=raster_Tmean_annual*rc
multiple_rainfall=multiple_rainfall*rc
raster_Tmin_annual=raster_Tmin_annual*rc
##########################################################################
load("data/Prosopis_glandulosa_presence_data")
##########################################################################

##########################################################################
predictors<-stack(raster_Tmean_annual,
                  raster_Tmin_annual,
                  raster_annualPPT,
                  multiple_rainfall
                  )

llmespg<-cbind(mesquite_pg_clean$lon[which(mesquite_pg_clean$year>1950)], mesquite_pg_clean$lat[which(mesquite_pg_clean$year>1950)])
llmespg<-llmespg[,-which(llmespg[,2]<0)]
presvals_pg<-as.data.frame(extract(predictors, llmespg))

names(presvals_pg)<-c("mean_annual_temperature","min_annual_temperature","total_annual_rainfall","multiple_rainfall")
names(predictors)<-names(presvals_pg)
##########################################################################

#training set (absences)
backgr <- dismo::randomPoints(predictors, 2000)
absvals <- as.data.frame(extract(predictors, backgr))
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
m1 <- glm(pg ~ mean_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m1_predict=predict(object = m1,newdata=as.data.frame(predictors))

m2 <- glm(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m2_predict=predict(object = m2,newdata=as.data.frame(predictors))


m3<-gam::gam(pg ~ mean_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m3_predict<-gam::predict.Gam(object = m3,newdata=as.data.frame(predictors))

m4<-gam::gam(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m4_predict<-gam::predict.Gam(object = m4,newdata=as.data.frame(predictors))
##Test to see the colinearity in the data
alias(m1)
#summary of results
summary(m1)
summary(m2)
summary(m3)
summary(m4)
##########################################################################



#evaluation
nn<-raster(x = predictors)
dist_m1<-setValues(x = nn,values = m1_predict)
dist_m2<-setValues(x = nn,values = m2_predict)
dist_m3<-setValues(x = nn,values = m3_predict)
dist_m4<-setValues(x = nn,values = m4_predict)

dist_m1=mask(dist_m1,USA_map)
dist_m2=mask(dist_m2,USA_map)
dist_m3=mask(dist_m3,USA_map)
dist_m4=mask(dist_m4,USA_map)

eval_m1<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m1,x=predictors)
eval_m2<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m2,x=predictors)
eval_m3<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m3,x=predictors)
eval_m4<-dismo::evaluate(p=llmespg[-traning,],a=backgr,model = m4,x=predictors)

##############################################################################################################
USA_map_full<-maptools::readShapeSpatial("C:/Users/abaezaca/Dropbox (ASU)/Design Principles in CIS/maps_USA/tl_2017_us_state/tl_2017_us_state")
top = 49.3457868 # north lat
left = -124.7844079 # west long
right = -66.9513812 # east long
bottom =  24.7433195 # south lat
USA_map<-raster::crop(USA_map_full,raster::extent(left,right,bottom,top))
USA_map@data$NAME<-factor(USA_map@data$NAME)
##############################################################################################################

png(filename = "output_maps/comparison_withTmin.png",width = 20,height = 20,res=300,units = "cm")
par(mfrow=c(2,2))
plot(dist_m1>eval_m1@t[which.max(eval_m1@TPR+eval_m1@TNR)],col=c("grey","red"),main="GLM without Tmin")
plot(USA_map,add=T)
plot(dist_m2>eval_m2@t[which.max(eval_m2@TPR+eval_m2@TNR)],col=c("grey","red"),main="GLM with Tmin")
plot(USA_map,add=T)
plot(dist_m3>eval_m3@t[which.max(eval_m3@TPR+eval_m3@TNR)],col=c("grey","red"),main="GAM without Tmin")
plot(USA_map,add=T)
plot(dist_m4>eval_m4@t[which.max(eval_m4@TPR+eval_m4@TNR)],col=c("grey","red"),main="GAM with Tmin")
plot(USA_map,add=T)

dev.off()