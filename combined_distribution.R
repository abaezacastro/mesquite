combined_distribution <-function(){
maps_present_GLM<-estimation_present(estimation_method = "GLM")

maps_past_GLM<-estimation_past(estimation_method = "GLM")

maps_present_GAM<-estimation_present(estimation_method = "GAM")

maps_past_GAM<-estimation_past(estimation_method = "GAM")


maps_present_Bioclim<-estimation_present(estimation_method = "BIOCLIM")

maps_past_Bioclim<-estimation_past(estimation_method = "BIOCLIM")


P50=(maps_past_GLM[[1]]> 0.5)*(maps_past_GAM[[1]]> 0.5)*(maps_past_Bioclim[[1]]> 0.5)
P2000=(maps_past_GLM[[2]]> 0.5)*(maps_past_GAM[[2]]> 0.5)*(maps_past_Bioclim[[2]]> 0.5)
P2100=(maps_past_GLM[[3]]> 0.5)*(maps_past_GAM[[3]]> 0.5)*(maps_past_Bioclim[[3]]> 0.5)

R2000=(maps_present_GLM[[1]]>0.5) * (maps_present_GAM[[1]]> 0.5)*(maps_present_Bioclim[[1]]> 0.5)
R2100=(maps_present_GLM[[2]]>0.5) * (maps_present_GAM[[2]] > 0.5) *(maps_present_Bioclim[[2]]>0.5)

P50=(P50==1)
P2000=(P2000==1)
P2100=(P2100==1)

R2000=(R2000==1)
R2100=(R2100==1)

maps_P<-stack(P50,P2000,P2100)
maps_R<- stack(R2000,R2100)
mapas<-list(maps_P,maps_R)
return(mapas)
}