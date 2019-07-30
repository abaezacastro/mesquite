load("raster_Tmean_current")
load("annual_total_rainfall_raster_current")
load("raster_Tmin_current")
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
load("Prosopis_glandulosa_presence_data")
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
backgr <- randomPoints(predictors, 2000)
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
#image(cor(sdmdata,use="pairwise.complete.obs"))
m1 <- glm(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m2 <- glm(pg ~ mean_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m3<-gam::gam(pg ~ mean_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
m4<-gam::gam(pg ~ mean_annual_temperature+min_annual_temperature+multiple_rainfall+total_annual_rainfall, data=traning_point,na.action = na.omit)
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
nn<-setValues(x = nn,values = p_training)
eval<-evaluate(p=llmespg[-traning,],a=backgr,model = m4,x=predictors)
data(wrld_simpl)




