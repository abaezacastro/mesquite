
estimation_present<-function(estimation_method ="GLM"){
metodo <- estimation_method

##########################################################################
#load("data/Prosopis_glandulosa_presence_data")
prosopis_species_data <- readRDS("data/species_presence_data")
##########################################################################

##########################################################################
#get cliamte scenario
predictors_current=get_current_climate()
##########################################################################

#llmespg<-cbind(mesquite_pg_clean$lon[which(mesquite_pg_clean$year>1950)], mesquite_pg_clean$lat[which(mesquite_pg_clean$year>1950)])
llmespg<-prosopis_species_data[which(prosopis_species_data$year>1950),1:2]
#llmespg<-llmespg[,-which(llmespg[,2]<0)]
presvals_pg<-as.data.frame(raster::extract(predictors_current, llmespg))
names(presvals_pg) <- names(predictors_current)

#backgr <- dismo::randomPoints(predictors_current, 2000) #if new background values are needed
#saveRDS(backgr,file = "data/background_points")
backgr<-readRDS("data/background_points")
bg_sam<-sample(x = 1:2000,size = 500,replace = F)
absvals <- as.data.frame(raster::extract(predictors_current, backgr[bg_sam,]))
names(absvals)<-names(presvals_pg)

pg <- c(rep(1, nrow(presvals_pg)), rep(0, nrow(absvals)))
sdmdata_pg <- data.frame(cbind(pg, rbind(presvals_pg, absvals)))

##########################################################################
valid_points=sdmdata_pg[which(sdmdata_pg$mean_annual_temperature!="NA"),]
traning<-sample(x = (1:length(valid_points$mean_annual_temperature)),size =round(length(valid_points$mean_annual_temperature)*0.75),replace = F)

traning_point=valid_points[traning,]

##########################################################################

#Estimations and predictions

if (metodo=="GLM"){

model <- glm(pg ~ mean_annual_temperature + min_annual_temperature+total_annual_rainfall + multiple_rainfall, data=traning_point,na.action = na.omit)
model_predict=predict(object = model,newdata=as.data.frame(predictors_current))
##########################################################################
#save parameters of the model
coefficients=coef(model)
summary_model=summary(model)$coefficients
write.csv(x=coefficients,file=paste(metodo,"_current_coefficients.csv"))
write.csv(x=summary_model,file=paste(metodo,"_current_summarymodel.csv"))
}
if (metodo=="GAM"){

model<-gam::gam(pg ~ mean_annual_temperature + min_annual_temperature+total_annual_rainfall + multiple_rainfall, data=traning_point,na.action = na.omit)
model_predict<-gam::predict.Gam(object = model,newdata=as.data.frame(predictors_current))

}


###########################################################
#bioclim uses only- presence data
if (metodo=="BIOCLIM"){

  valid_points=presvals_pg[which(presvals_pg$mean_annual_temperature!="NA"),]
  traning<-sample(x = (1:length(valid_points$mean_annual_temperature)),size =round(length(valid_points$mean_annual_temperature)*0.75),replace = F)
  traning_point=valid_points[traning,]
  
  model<-dismo::bioclim(traning_point[,-3])
  model_predict<-predict(model,predictors_current)

}

##########################################################################
#create raster with model prediction
if (metodo!="BIOCLIM"){
  nn<-raster(x = predictors_current)
  dist<-setValues(x = nn,values = model_predict)
}
else{
  dist <- model_predict
}
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
  eval <- dismo::evaluate(p=llmespg[-traning,],a=backgr[bg_sam,],model = model,x=predictors_current)

  predictors_future <- get_future_climate()
  nn_future <- raster(x = predictors_future)

  P2000 <- dist>eval@t[which.max(eval@TPR+eval@TNR)]

  if(metodo!="BIOCLIM"){
    model_predict_future <- predict(object = model,newdata=as.data.frame(predictors_future))
    dist_future <- setValues(x = nn_future,values = model_predict_future)
    
    }
  else{
    model_predict_future <- predict(model,predictors_future)
    dist_future <-model_predict_future
      }
  
  dist_future <-mask(dist_future,USA_map)
  P2100 <- dist_future>eval@t[which.max(eval@TPR+eval@TNR)]

  desagregation_factor=res(predictors_future)/res(predictors_current)  
  P2100 %>% disaggregate(fact=desagregation_factor) -> P2100
  crs(P2100) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84"
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

