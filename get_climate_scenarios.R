#############################################################################################
get_past_climate <- function(){
  load("data/raster_Tmean_historic")
  load("data/annual_total_rainfall_raster_historico")
  load("data/raster_Tmin_historic")
  multiple_rainfall=PPT_m_H*PPT_m_H
  ##########################################################################
  #re-clasification
  m <- c(-1, 0, NA,  0.1, 5000, 1)
  rclmat <- matrix(m, ncol=3, byrow=TRUE)
  # rm(m)
  
  rc <- reclassify(PPT_m_H, rclmat)
  
  
  total_annual_rainfall=PPT_m_H*rc
  raster_Tmean_annual_H=raster_Tmean_annual_H*rc
  multiple_rainfall=multiple_rainfall*rc
  raster_Tmin_annual_H=raster_Tmin_annual_H*rc
  
  predictors_past<-stack(raster_Tmean_annual_H,
                         raster_Tmin_annual_H,
                         total_annual_rainfall,
                         multiple_rainfall
  )
  
  names(predictors_past)<-c("mean_annual_temperature",
                            "min_annual_temperature",
                            "total_annual_rainfall",
                            "multiple_rainfall")
  
  return(predictors_past)
 
}

#############################################################################################

get_current_climate <- function(){
  load("data/raster_Tmean_current")
  load("data/annual_total_rainfall_raster_current")
  load("data/raster_Tmin_current")
  raster_annualPPT<-calc(raster_annualPPT,fun=function(x){x/10})
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
  
  predictors_current<-stack(raster_Tmean_annual,
                            raster_Tmin_annual,
                            total_annual_rainfall,
                            multiple_rainfall
  )
  
  
  names(predictors_current)<-c("mean_annual_temperature",
                               "min_annual_temperature",
                               "total_annual_rainfall",
                               "multiple_rainfall")
  return(predictors_current)
}

#############################################################################################
#future climate

get_future_climate <- function(save_file =FALSE){
  R_fp_r <- stack("c:/Users/abaezaca/Dropbox (ASU)/Mesquite/Modeling2/future_scenarios/Stats_pr.nc",varname="pr_pm")
  
  T_fp_r <- stack("c:/Users/abaezaca/Dropbox (ASU)/Mesquite/Modeling2/future_scenarios/Stats_tas.nc",varname="tas_pm")
  
  
  T_min<-raster("c:/Users/abaezaca/Dropbox (ASU)/Mesquite/Modeling2/future_scenarios/Stats_tasmin.nc",varname="tasmin_pm")
  
  
  #Change Projection
  
  load("c:/Users/abaezaca/Documents/GitHub/mesquite/data/raster_Tmean_current")
  
  
  T_fp_r %>% rotate() %>% crop(y=raster_Tmean_annual)-> Tmean
  R_fp_r %>% rotate() %>% crop(y=raster_Tmean_annual) -> rain
  T_min %>% rotate() %>% crop(y=raster_Tmean_annual) -> Tmin
  
  predictors_future<-stack(Tmean,Tmin,rain,rain*rain)
  names(predictors_future)=c("mean_annual_temperature","min_annual_temperature","total_annual_rainfall","multiple_rainfall")
  
  if(save_file == TRUE){
  saveRDS(predictors_future,file = "climate_scenario_Tmin_Tmean_rain_2080-2100")
  }
  return(predictors_future)
}
#############################################################################################



