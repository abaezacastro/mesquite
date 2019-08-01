
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
saveRDS(predictors_future,file = "climate_scenario_Tmin_Tmean_rain_2080-2100")
rm(Tmean)
rm(Tmin)
rm(rain)
rm(T_fp_r)
rm(R_fp_r)
rm(T_min)
rm(raster_Tmean_annual)




