
get_future_climate <- function(save_file =FALSE){
R_fp_r <- stack("data/futureScenarios/Stats_pr.nc",varname="pr_pm")
                
T_fp_r <- stack("data/futureScenarios/Stats_tas.nc",varname="tas_pm")


T_min<-raster("data/futureScenarios/Stats_tasmin.nc",varname="tasmin_pm")


#Change Projection

load("data/raster_Tmean_current")


T_fp_r %>% rotate() %>% crop(y=raster_Tmean_annual)-> Tmean
R_fp_r %>% rotate() %>% crop(y=raster_Tmean_annual) -> rain
T_min %>% rotate() %>% crop(y=raster_Tmean_annual) -> Tmin

predictors_future<-stack(Tmean,Tmin,rain,rain*rain)
names(predictors_future)=c("mean_annual_temperature","min_annual_temperature","total_annual_rainfall","multiple_rainfall")
saveRDS(predictors_future,file = "climate_scenario_Tmin_Tmean_rain_2080-2100")
rm(Tmean)
rm(Tmin)
rm(rain)
rm(T_fp_r)
rm(R_fp_r)
rm(T_min)
rm(raster_Tmean_annual)
return(predictors_future)
}


