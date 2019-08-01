
estimation_present<-function(estimation_method ="GLM"){
metodo <- estimation_method

##########################################################################
load("data/Prosopis_glandulosa_presence_data")
##########################################################################

##########################################################################
#get cliamte scenario
predictors_current=get_current_climate()
##########################################################################

llmespg<-cbind(mesquite_pg_clean$lon[which(mesquite_pg_clean$year>1950)], mesquite_pg_clean$lat[which(mesquite_pg_clean$year>1950)])
llmespg<-llmespg[,-which(llmespg[,2]<0)]
presvals_pg<-as.data.frame(raster::extract(predictors_current, llmespg))
names(presvals_pg)<-c("mean_annual_temperature","min_annual_temperature","total_annual_rainfall","multiple_rainfall")

#backgr <- dismo::randomPoints(predictors_current, 2000) #if new background values are needed
#saveRDS(backgr,file = "data/background_points")
backgr<-readRDS("data/background_points")
absvals <- as.data.frame(raster::extract(predictors_current, backgr))
names(absvals)<-names(presvals_pg)

pg <- c(rep(1, nrow(presvals_pg)), rep(0, nrow(absvals)))
sdmdata_pg <- data.frame(cbind(pg, rbind(presvals_pg, absvals)))

##########################################################################
valid_points=which(sdmdata_pg$mean_annual_temperature!="NaN")
traning<-sample(x = (1:length(sdmdata_pg$mean_annual_temperature))[valid_points],size =round(length(valid_points)/2),replace = F)

traning_point=sdmdata_pg[traning,]

##########################################################################

#Estimations and predictions

if (metodo=="GLM"){

model <- glm(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
model_predict=predict(object = model,newdata=as.data.frame(predictors_current))

}
if (metodo=="GAM"){

model<-gam::gam(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
model_predict<-gam::predict.Gam(object = m4,newdata=as.data.frame(predictors_current))

}


###########################################################
#bioclim uses only- presence data
if (metodo=="BIOCLIM"){

valid_points=which(presvals_pg$mean_annual_temperature!="NaN")
traning<-sample(x = (1:length(presvals_pg$mean_annual_temperature))[valid_points],size =round(length(valid_points)/2),replace = F)
traning_point=presvals_pg[traning,]

model<-dismo::bioclim(traning_point)
model_predict<-predict(m6,predictors_current)

}


##########################################################################
#create raster with model prediction
nn<-raster(x = predictors_current)
dist<-setValues(x = nn,values = model_predict)

##############################################################################################################
USA_map_full<-maptools::readShapeSpatial("C:/Users/abaezaca/Dropbox (ASU)/Design Principles in CIS/maps_USA/tl_2017_us_state/tl_2017_us_state")
top = 49.3457868 # north lat
left = -124.7844079 # west long
right = -66.9513812 # east long
bottom =  24.7433195 # south lat
USA_map<-crop(USA_map_full,raster::extent(left,right,bottom,top))
USA_map@data$NAME<-factor(USA_map@data$NAME)
rm(top)
rm(left)
rm(right)
rm(bottom)
##############################################################################################################
#evaluation
  dist <- mask(dist,USA_map)
  eval <- dismo::evaluate(p=llmespg[-traning,],a=backgr,model = model,x=predictors_current)

  predictors_future <- get_future_climate()
nn_future <- raster(x = predictors_future)

  P2000 <- dist>eval@t[which.max(eval@TPR+eval@TNR)]
  model_predict_future <- predict(object = model,newdata=as.data.frame(predictors_future))
  dist_future <- setValues(x = nn_future,values = model_predict_future)
  dist_future <-mask(dist_future,USA_map)
  P2100 <- dist_future>eval@t[which.max(eval@TPR+eval@TNR)]

  P2100%>%
    disaggregate(fact=res(predictors_future)/res(predictors_current)) -> P2100
  crs(P2100) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
    
  P2000 %>% crop(P2100) -> P2000

  
  crs(P2000) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
  
  P2100 <- resample(P2100,P2000)

return(stack(P2000,P2100))  
rm(dist)
rm(dist_future)
rm(nn_future)
rm(sdmdata_pg)
rm(presvals_pg)
rm(valid_points)
rm(llmespg)
rm(nn)
rm(eval)
rm(model_predict_future)
rm(backgr)
rm(traning_point)
rm(model_predict)
rm(absvals)
rm(mesquite_pg_clean)
rm(model)
rm(rclmat)
rm(USA_map_full)
rm(pg)
rm(traning)

}
###############################################################################3
# png(filename = "output_maps/comparison_withTmin.png",width = 20,height = 30,res=300,units = "cm")
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



