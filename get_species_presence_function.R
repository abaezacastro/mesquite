get_species_presence <-function(especie='Prosopis glandulosa'){
  
  mesquite_p<-dismo::gbif(genus = 'Prosopis',species = especie,geo = TRUE)

  #mesquite_p_NNA<-mesquite_p[which(is.na(mesquite_p$lat)==F & is.na(mesquite_p$lon)==F),]
  #remove duplicated data
  #mesquite_p_clean<- mesquite_p_NNA[which(duplicated(mesquite_p_NNA[, 1:10])==FALSE),]
  return(mesquite_p)
}
load("data/Prosopis_glandulosa_presence_data")
prosopis_species_data_pg <- mesquite_pg_clean[,c("lon","lat","year")]
prosopis_species_data_pj <- get_species_presence('Prosopis juliflora')
prosopis_species_data_pv <- get_species_presence('Prosopis velutina')


prosopis_species_data <- rbind(
  prosopis_species_data_pg[which(is.na(prosopis_species_data_pg$lat)==F & is.na(prosopis_species_data_pg$lon)==F),c("lon","lat","year")],
  prosopis_species_data_pj[which(is.na(prosopis_species_data_pj$lat)==F & is.na(prosopis_species_data_pj$lon)==F),c("lon","lat","year")],
  prosopis_species_data_pv[which(is.na(prosopis_species_data_pv$lat)==F & is.na(prosopis_species_data_pv$lon)==F),c("lon","lat","year")]
)

saveRDS(prosopis_species_data,file = "data/species_presence_data")
