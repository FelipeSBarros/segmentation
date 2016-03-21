segmentation <- function(
  envLayer = predictors, #raster Layer or raster stack
  studyArea = studyArea, # SpatialPolygonsDataFrame
  projName = 'MT', # Sufix to be used when saving results
  randomforest = TRUE,
  random.pt = NULL, # Number of random points to be genarated to run randomForest
  Kmeans = TRUE,
  ngroup = NULL,
  polygonize = FALSE
  
){
  library(rgdal)
  library(raster)
  library(rgeos)
  library(randomForest)
  library(ggplot2)
  if (randomforest){
    #Generating random points:
    if (is.null(random.pt)){
      random.pt <- as.integer(0.01 * ncell(envLayer)) 
    }
    cat("Generating", random.pt,"random points")
    set.seed(2602)
    sp.pts <- spsample(studyArea, random.pt, type="random")
    
    #Extracting raster values for each point
    sp.pts$prec <- extract(envLayer, sp.pts)
    
    #removing NA
    if (any(is.na(sp.pts@data))){
      sp.pts <- sp.pts[-which(is.na(sp.pts@data)),]
      }
    
    #creating object with the raster value to run kmeans
    mydata <- as.data.frame(sp.pts@data)
    
    #mydata <- scale(mydata) # standardize variables 
    
    # Determine number of clusters
    wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
    
    for (i in 2:(as.integer(0.01 * random.pt))) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
    wss2 <- as.data.frame(wss)
    wss2$row <- as.integer(rownames(wss2))
    
    # Saving plot
    ggplot(wss2, aes(x = row, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size=17))
    ggsave(paste0("Kmeans_clusterAnalysis_",projName, ".png"), dpi=300)
    dev.off()
    
    plot(wss2$row, wss2$wss, xlab = "Number of Clusters", ylab = "Within groups sum of squares", type='b')
    
    if (is.null(ngroup)){
      # Number of classes to be classified
      ngroup <- readline("How many groups should be created?")  
    }
    
    # K-Means Cluster Analysis
    fit <- kmeans(mydata, ngroup, algorithm="Lloyd")
    
    # append cluster assignment
    mydata <- data.frame(mydata, fit$cluster) 
    
    # CREATE RF MODEL
    names(envLayer) <- paste0(rep('band', nlayers(envLayer)),1:nlayers(envLayer))
    names(mydata) <- c(paste0(rep('band', nlayers(envLayer)),1:nlayers(envLayer)),'fit.cluster')
    
    rf.mdl <- randomForest (as.factor(mydata$fit.cluster)~.,data=mydata)
    
    rf_segmentation <- raster::predict(envLayer, rf.mdl, progress="text", type='response')
    
    if(polygonize){
      cat("converting and saving to Polygons")
      # Saving Vector output (Polygon)
      rf_segmentation <- rasterToPolygons(rf_segmentation, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
      writeOGR(rf_segmentation, dsn='./', layer=paste0('rf_segmentation_',projName), driver ="ESRI Shapefile", overwrite_layer=TRUE)
    } else {
      cat("Saving raster result")
      #Saving RASTER output
      writeRaster(rf_segmentation, filename=paste0('./rf_segmentation_',projName,'.tif'), overwrite=TRUE)    }
  }
  if (Kmeans){
    # Creating mask to remove background
    mask <- envLayer[[1]]>=0
    
    # Removing NA from background
    envLayer[is.na(envLayer)] <- 999
    
    # Creating data to Kmeans analysis
    mydata <- values(envLayer)
    
    # mydata <- scale(mydata) # standardize variables 
    
    # Determine number of clusters
    #wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))
    
    #for (i in 2:(as.integer(0.00001 * ncell(envLayer)))) wss[i] <- sum(kmeans(mydata, centers=i)$withinss)
    #wss2 <- as.data.frame(wss)
    #wss2$row <- as.integer(rownames(wss2))
    
    # Saving plot
    #ggplot(wss2, aes(x = row, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size=17))
    #ggsave(paste0("Kmeans_clusterAnalysis_",projName, ".png"), dpi=300)
    
    #ggplot(wss2, aes(x = row, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size=17))
    if (is.null(ngroup)){
      # Number of classes to be classified
      ngroup <- readline("How many groups should be created?")  
    }

    # K-Means Cluster Analysis
    fit <- kmeans(mydata, (as.integer(ngroup) + 1), algorithm="Lloyd") # n cluster solution
    
    # append cluster assignment
    mydata <- data.frame(mydata, fit$cluster) 
    
    # creating a raster layer to recieve goup values
    km_SegmentationRaster <- envLayer[[1]]
    # changing pixel values to group values
    values(km_SegmentationRaster) <- mydata[,ncol(mydata)]
    
    #making to remove background
    km_SegmentationRaster <- mask(km_SegmentationRaster, mask)
    
    if(polygonize){
      # Saving Vector output (Polygon)
      km_SegmentationRaster <- rasterToPolygons(km_SegmentationRaster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE)
      writeOGR(km_SegmentationRaster, dsn='./', layer=paste0('km_segmentation_',projName), driver ="ESRI Shapefile", overwrite_layer=TRUE)
    } else {
      #Saving RASTER output
      writeRaster(km_SegmentationRaster, filename=paste0('./km_segmentation_',projName,'.tif'), overwrite=TRUE) }
  }
  }