satPrep <- function(rasterFolder = './ModData/1985/',
                    pattern='.TIF$',
                    crop.Ext = NULL,
                    ndvi = TRUE,
                    evi = TRUE,
                    savi = TRUE,
                    nameOutput = 'stack_1985'
){
  library(raster)
  #Stack of layers
  satImagery <- stack(list.files( rasterFolder, pattern=pattern, full.names = TRUE))
  
  if (!is.null(crop.Ext)){
    if (! proj4string(crop.Ext) == proj4string(satImagery)){
      crop.Ext <- spTransform(crop.Ext, CRS(proj4string(satImagery)))
    }
    
    satImagery <- crop(satImagery, crop.Ext)
  }
  
  if (ndvi){
    # NDVI----
    ndvi <- overlay(satImagery[[4]], satImagery[[3]], fun=function(x,y){(x-y)/(x+y)})
    names(ndvi) <- 'ndvi'
    satImagery <- addLayer(satImagery, ndvi)
  }
  
  if (evi){
    # EVI ----
    evi <- overlay(satImagery[[4]], satImagery[[3]], satImagery[[1]], fun=function(x,y,z){(2.5*(x-y)/((x+6)*(y-7.5)*(z+1)))})
    names(evi) <- 'evi'
    satImagery <- addLayer(satImagery, evi)
  }
  
  if(savi){
    # SAVI ----
    savi <- overlay(satImagery[[5]], satImagery[[4]], fun=function(x,y){((x-y)/(x+y+0.5))*(1.5)})
    names(savi) <- 'savi'
    satImagery <- addLayer(satImagery, savi)
  }
  writeRaster(satImagery, filename = paste0(rasterFolder,nameOutput,'.tif'), overwrite = TRUE)
  return(satImagery)
}

ForIndex <- function(
  studyArea = readOGR(dsn = './shp/', layer='DoisIrmaos'),
  idField = 'SETOR',
  group = 'OBRA',
  satFolder = './ModData/1985/',
  rstPattern = 'stack',
  idRef = 12,
  stats='mean',
  rastLayer=c(8, 9, 10),
  respRatio=TRUE
){
  #This function is based only on spectral distance or mean. It is not needed a R output from Segmentation_fct
  
  library(raster)
  library(rgdal)
  library(rgeos)
  
  studyArea@data$AREA <- gArea(studyArea,byid=T)
  studyArea@data$id_fct <- 1:nrow(studyArea)
  if (is.null(group)) {
    styduAreaDF <- as.data.frame(studyArea[,c(idField, 'id_fct', 'AREA')])
    } else {
    styduAreaDF <- as.data.frame(studyArea[,c(idField, 'id_fct', group, 'AREA')])  
  }
  
  # Ordering by idField #Not necessari anymore
  #styduAreaDF <- styduAreaDF[order(styduAreaDF[, idField]),]
  
  # applying rasterStack function to raster 1985 ----
  sat <- stack(list.files(paste0(satFolder), pattern=paste0(rstPattern, '.*tif$'), full.names = TRUE))
  #sat <- sat[[1]]
  
  # projecting sutydArea if necessary
  if (!identicalCRS(sat, studyArea)){
    # Alterando CRS: studyArea
    studyArea <- spTransform(studyArea, sat@crs)
    }
  
  # Rasterizing studyArea
  rasterArea <- rasterize(studyArea, sat, field=studyArea@data[,'id_fct'])
  # plot(rasterArea, ext=extent(studyArea))
  
  # Estimating the 'stast' defined by parameter of layers reference (NDVI, SAVI, ...) for each idField 'zone'
  for (a in rastLayer){
    # Estimando valor a partir da estatistica definida por parametro
    rastStats <- as.data.frame(zonal(sat[[a]], rasterArea, fun=stats))
    rastStats <- merge(styduAreaDF, rastStats, by.x='id_fct', by.y='zone', all = TRUE)[ncol(styduAreaDF)+1]
    styduAreaDF[, paste0(stats, a)] <- rastStats
    
    #Standart deviation
    rastStats2 <- as.data.frame(zonal(sat[[a]], rasterArea, fun=sd))
    rastStats2 <- merge(styduAreaDF, rastStats2, by.x='id_fct', by.y='zone', all = TRUE)[length(styduAreaDF)+1]
    styduAreaDF[, paste0('sd', a)] <- rastStats2
    
    #If should use response ratio, must have idRef defined!
    if (respRatio & !is.null(idRef)){
      #Changing idref provided by parameters for the one created in function
      idRef <- styduAreaDF[which(styduAreaDF$iis_id %in% idRef),'id_fct']
      
      styduAreaDF[ , ((ncol(styduAreaDF))+1) ] <- 1
      
      for (b in (styduAreaDF[#-c(idRef)
        , 'id_fct'])){
        styduAreaDF[b, ncol(styduAreaDF) ] <- log(rastStats[b, stats]/rastStats[idRef, stats])
        # Estimando distancia entre areas de estudo e a area de ref
      }
      names(styduAreaDF)[ncol(styduAreaDF)] <- paste0("RR", a)
    }
    # In case respRatio is false and idRef provided, run dist from study area and ref area (idRef)
    if (!respRatio & !is.null(idRef)){
      styduAreaDF[ , ((ncol(styduAreaDF))+1) ] <- 1
      
      idRef <- styduAreaDF[which(styduAreaDF$iis_id %in% idRef),'id_fct']
      dist <- as.matrix(dist(rastStats[,stats], method = "euclidean"))
      
      for (b in (styduAreaDF[, 'id_fct'])){
        #styduAreaDF[b, paste0('mean',a)] <- rastStats[b, stats]
        #dist[rastStats[b, 'zone'], idRef]
        styduAreaDF[b, ncol(styduAreaDF)] <- dist[b, idRef]
      }
      names(styduAreaDF)[ncol(styduAreaDF)] <- paste0("Dist", a)
      
    }
    # Removendo polgonos menores que 4 pixels
    
  }
  styduAreaDF$data <- satFolder
  
  #styduAreaDF <- styduAreaDF[-which(styduAreaDF$AREA<3600),] #4 pixels
  
  return(styduAreaDF)
}

ForIndexZone <- function(
  studyArea = readOGR(dsn = './shp/', layer='DoisIrmaos'),
  idField = 'SETOR',
  group = 'OBRA',
  satFolder = './ModData/1985/',
  idRef = NULL,
  stats='mean'
){
  #This function is based on the unsupervised results. Must have a R output from Segmentation_fct
  library(raster)
  library(rgdal)
  library(rgeos)
  library(RColorBrewer)
  pallete <- brewer.pal(8, 'Accent')
  studyArea@data$AREA <- gArea(studyArea,byid=T)
  if (is.null(group)) styduAreaDF <- as.data.frame(studyArea[,c(idField, 'AREA')]) else {
    styduAreaDF <- as.data.frame(studyArea[,c(idField, group, 'AREA')])  
  }
  
  # Ordering by idField
  styduAreaDF <- styduAreaDF[order(styduAreaDF[, idField]),]
  
  # applying rasterStack function to raster 1985 ----
  sat <- stack(list.files(paste0(satFolder), pattern='.tif$', full.names = TRUE))
  
  if (!identicalCRS(sat, studyArea)){
    # Alterando CRS: studyArea
    studyArea <- spTransform(studyArea, sat@crs)
  }
  
  # Rasterizando Area de estudo
  rasterArea <- rasterize(studyArea, sat, field=studyArea@data$SETOR)
  #plot(rasterArea, col=pallete)
  
  #Carregando classificação
  rastClass <- raster(list.files(path = paste0(satFolder), pattern='fuzzy_segmentation', full.names = TRUE))
  
  plot(rastClass, col=pallete)
  idRef <- readline("Which group must be used as reference?")
  
  # Estimando valor mais frequente
  #rastStats <- zonal(sat[[rastLayer]], rasterArea, fun=stats)
  zonalClass <- zonal(rastClass, rasterArea, fun='modal')
  
  #Carregando estatisticas
  load(list.files(path = paste0(satFolder), pattern='_Fuzzy_.RData', full.names = TRUE))
  
  #Estimando distancia entre areas de estudo
  dist <- as.matrix(dist(fuzzy$centers, method = "euclidean"))
 
  #Criando data frame 
  for (a in 1:nrow(styduAreaDF)){
    styduAreaDF[a, 'modal'] <- zonalClass[a, 'modal']
    styduAreaDF[a, 'dist'] <- dist[zonalClass[a, 'modal'], idRef]
  }
  
  return(styduAreaDF) 
}