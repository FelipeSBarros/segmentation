segmentation <- function(envLayer = envLayer, #raster Layer or raster stack
                         studyArea = studyArea, # SpatialPolygonsDataFrame
                         projName = 'MT', # Sufix to be used when saving results
                         randomforest = TRUE,
                         Kmeans = TRUE,
                         fuzzy.cluster = FALSE,
                         random.pt = NULL, # Number of random points to be genarated to run randomForest
                         ngroup = NULL,
                         save.shp = FALSE,
                         save.raster = FALSE,
                         save.plot=FALSE,
                         seed = 123) {
  library(rgdal)
  library(raster)
  library(rgeos)
  library(randomForest)
  library(ggplot2)

  if (randomforest) {
    #Generating random points:
    if (is.null(random.pt)) {
      random.pt <- as.integer(0.01 * ncell(envLayer))
    }
    cat("Generating", random.pt,"random points \n")
    set.seed(seed)
    sp.pts <- spsample(studyArea, random.pt, type = "random")
    
    #Extracting raster values for each point
    sp.pts$prec <- extract(envLayer, sp.pts)
    
    #removing NA
    if (any(is.na(sp.pts@data))) {
      sp.pts <- sp.pts[-which(is.na(sp.pts@data)),]
    }
    
    #creating object with the raster value to run kmeans
    mydata <- as.data.frame(sp.pts@data)
    
    # Not used for this prupouse
    #mydata <- scale(mydata) # standardize variables
    
    # Determine number of clusters
    wss <- (nrow(mydata) - 1) * sum(apply(mydata,2,var))
    
    for (i in 2:(as.integer(0.01 * random.pt)))
      wss[i] <- sum(kmeans(mydata, centers = i)$withinss)
    wss2 <- as.data.frame(wss)
    wss2$row <- as.integer(rownames(wss2))
    
    # Saving plot
    ggplot(wss2, aes(x = row, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size = 17))
    ggsave(paste0("Kmeans_clusterAnalysis_",projName, ".png"), dpi = 300)
    dev.off()
    
    plot(wss2$row, wss2$wss, xlab = "Number of Clusters", ylab = "Within groups sum of squares", type = 'b')
    
    if (is.null(ngroup)) {
      # Number of classes to be classified
      ngroup <- readline("How many groups should be created?")
    }
    # K-Means Cluster Analysis
    fit <- kmeans(
      mydata, ngroup, algorithm = c("Hartigan-Wong",
                                    "Lloyd",
                                    "Forgy",
                                    "MacQueen")[1])
    
    # append cluster assignment
    mydata <- data.frame(mydata, fit$cluster)
    
    # Creating randomForest Model
    names(envLayer) <-
      paste0(rep('band', nlayers(envLayer)),1:nlayers(envLayer))
    names(mydata) <-
      c(paste0(rep('band', nlayers(envLayer)),1:nlayers(envLayer)),'fit.cluster')
    
    rf.mdl <-
      randomForest (as.factor(mydata$fit.cluster) ~ .,data = mydata)
    
    rf_segmentation <-
      raster::predict(envLayer, rf.mdl, progress = "text", type = 'response')
    
    if (polygonize) {
      cat("converting and saving to Polygons \n")
      # Saving Vector output (Polygon)
      rf_segmentation <-
        rasterToPolygons(
          rf_segmentation, fun = NULL, n = 4, na.rm = TRUE, digits = 12, dissolve = FALSE)
      writeOGR(
        rf_segmentation, dsn = './', layer = paste0('rf_segmentation_',projName), driver =
          "ESRI Shapefile", overwrite_layer = TRUE
      )
    } else {
      cat("Saving raster result \n")
      #Saving RASTER output
      writeRaster(
        rf_segmentation, filename = paste0('./rf_segmentation_',projName,'.tif'), overwrite = TRUE
      )
    }
    cat("randomForest segmentation done. \n")
  }
  if (Kmeans) {
    mydata <- getValues(envLayer)
    
    # Identificando pixel nao NA
    a <- which(!is.na(mydata))
    
    # Creating data to Kmeans analysis
    mydata <- na.omit(values(envLayer))
    
    # Not used for this prupouse
    mydata <- scale(mydata) # standardize variables
    
    if(save.plot){  
      
      # Determine number of clusters
      wss <- (nrow(mydata) - 1) * sum(apply(mydata,2,var))
    
      for (b in 2:20)
      wss[b] <- sum(kmeans(mydata, centers = i)$withinss)
      wss2 <- as.data.frame(wss)
      wss2$row <- as.integer(rownames(wss2))
      
      #plotting results
      plot(wss2$row, wss2$wss, xlab = "Number of Clusters", ylab = "Within groups sum of squares", type = 'b')
      
      # Saving plot
      ggplot(wss2, aes(x = row, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size = 17))
      ggsave(paste0("Kmeans_clusterAnalysis_",projName, ".png"), dpi = 300)
      dev.off()
    }
    
    if (is.null(ngroup)) {
      # Number of classes to be classified
      ngroup <- readline("How many groups should be created?")
      }
      
    # K-Means Cluster Analysis
    set.seed(seed)
    fit <-
      kmeans(
        mydata, ngroup, algorithm = c("Hartigan-Wong",
                                      "Lloyd",
                                      "Forgy",
                                      "MacQueen")[1], nstart=10)
    
    # creating a raster layer to recieve goup values
    km_SegmentationRaster <- envLayer[[1]]
    
    # changing pixel values to group values
    km_SegmentationRaster[a] <- fit$cluster
    
    # plot(mydata[,c(7:7)], col=fit$cluster, pch='.')
    # points(fit$centers[,c(1, 2)], col=1:3, pch=8, cex=2)
    
    if (save.shp) {
        # Saving Vector output (Polygon)
        km_SegmentationRaster <-
          rasterToPolygons(
            km_SegmentationRaster, fun = NULL, n = 4, na.rm = TRUE, digits = 12, dissolve =  FALSE)
        writeOGR(
          km_SegmentationRaster, dsn = './', layer = paste0('km_segmentation_',projName), driver = "ESRI Shapefile", overwrite_layer = TRUE)
      } 
    if (save.raster) {
        #Saving RASTER output
        cat("Saving raster result \n")
        writeRaster(
          km_SegmentationRaster, filename = paste0('./km_segmentation_',projName,'.tif'), overwrite = TRUE)
        }
    cat("Kmeans segmentation done. \n")
  }
  if (fuzzy.cluster){
    library(e1071)
    
    mydata <- getValues(envLayer)
    
    # Identificando pixel nao NA
    a <- which(!is.na(mydata))
    
    # Creating data to Kmeans analysis
    mydata <- na.omit(values(envLayer))
    
    # Not used for this prupouse
    mydata <- scale(mydata) # standardize variables
    
    set.seed(seed)
    fuzzy <- cmeans(mydata, ngroup, iter.max = 100, verbose = FALSE,
                    dist = c("euclidean", "manhattan")[1], method = c("cmeans","ufcl")[1], m = 2,
                    rate.par = NULL, weights = 1, control = list())
    
    # creating a raster layer to recieve goup values
    fuzzy_SegmentationRaster <- envLayer[[1]]
    
    # changing pixel values to group values
    fuzzy_SegmentationRaster[a] <- fuzzy$cluster
    
    if (save.shp) {
      # Saving Vector output (Polygon)
      fuzzy_SegmentationRaster <-
        rasterToPolygons(
          fuzzy_SegmentationRaster, fun = NULL, n = 4, na.rm = TRUE, digits = 12, dissolve =  FALSE)
      writeOGR(
        fuzzy_SegmentationRaster, dsn = './', layer = paste0('fuzzy_segmentation_',projName), driver = "ESRI Shapefile", overwrite_layer = TRUE)
    } 
    if (save.raster) {
      #Saving RASTER output
      cat("Saving raster result \n")
      writeRaster(
        fuzzy_SegmentationRaster, filename = paste0('./fuzzy_segmentation_',projName,'.tif'), overwrite = TRUE)
    }
    cat("Fuzzy segmentation done. \n")
  }
  }