precAnalysis <- function(envLayer = MT_prec,
                zones = MT_20,
                stats = c('mean','sd'),
                projName = 'MT_yr_prec'){
                  
  ##Runing statistics for each zone:
  mean <- as.data.frame(zonal(envLayer, zones, fun=stats[1]))
  sd <- as.data.frame(zonal(envLayer, zones, fun=stats[2]))
  #head(mean)
  #head(sd)
  
  ##Saving the statistics results:
  if (!file.exists("./out")){
    dir.create("./out")}
  write.csv(mean, file = paste0('./out/', stats[1], '_', projName,'.csv'))
  write.csv(sd, file = paste0('./out/', stats[2], '_', projName,'.csv'))
  
  ## Computing anual precipitation
  mean$total <- rowSums(mean[,2:length(mean)], na.rm = TRUE)
  #sd$total <- rowSums(sd[,2:length(sd)], na.rm = TRUE)

  ##Reclassifying map to have anual mean precipitation
  mean_map <- subs(zones, mean, by='zone', which='total')
  plot(mean_map)
  writeRaster(mean_map, filename = paste0('./out/', stats[1], '_', projName, '.tif'), overwrite=TRUE)
}