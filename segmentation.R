segmentation <- function(envLayer = envLayer, #raster Layer or raster stack
                         studyArea = studyArea, # SpatialPolygonsDataFrame
                         projName = 'MT', # Sufix to be used when saving results
                         folder = './',
                         randomforest = FALSE,
                         Kmeans = FALSE,
                         fuzzy.cluster = FALSE,
                         random.pt = NULL, # Number of random points to be genarated to run randomForest
                         ngroup = NULL,
                         save.shp = FALSE,
                         save.raster = FALSE,
                         explore=TRUE,
                         h.life=TRUE,
                         save.fit=TRUE,
                         seed = 123) {
  library(rgdal)
  library(raster)
  library(rgeos)
  library(randomForest)
  library(ggplot2)

  if(explore){  
    
    # Will analyse the whithin sum of squares variation (cahnges) for different group number
    cat("Processing K-means for exploratory analysis  \n")
    mydata <- getValues(envLayer)
    
    # Identificando pixel nao NA
    a <- which(!is.na(mydata))
    
    # Creating data to Kmeans analysis
    mydata <- na.omit(values(envLayer))
    
    # Not used for this prupouse
    mydata <- scale(mydata) # standardize variables
    
    set.seed(seed)
    wss <- matrix(ncol=3, nrow=2*ngroup)
    wss[1,1] <- (nrow(mydata) - 1) * sum(apply(mydata,2,var))
    wss[1,2] <- 0
    wss[1,3] <- 0
    
    for (b in 2:nrow(wss)){
      set.seed(seed)
      wss[b, 1] <- sum(kmeans(mydata, centers = b)$withinss)
      wss[b, 2] <- abs(wss[b-1, 1] - wss[b, 1])
      wss[b, 3] <- wss[b, 2] / wss[2, 2]
      #wss[b, 4] <- wss[b, 1] / wss[1, 1]
    }
    
    wss <- as.data.frame(wss)
    names(wss) <- c('wss', 'DeltaChange', 'RatioChange')
    wss$group <- as.integer(rownames(wss))
    if (h.life){
      # will analyse the amount of groups needed to achieve the half exponential decay (half life)
      # exp. decay fit
      exp.decay = lm(log(wss$wss) ~ wss$group)
      
      # exp. decay function using the fitted parameters
      fun = function(x){exp(exp.decay$coefficients[2]*x + exp.decay$coefficients[1])}
      
      # exp. decay mean lifetime (scaling time)
      tau = abs(1/exp.decay$coefficients[2])
      
      # exp. decay half life
      # quantos grupos são necessários para reducao à metade do variacao inicial.
      h.life = tau * log(2)
      # log(2)/teste$coefficients[2] #same as previous
      
      # visualising data and fitted function
      ggplot(wss, aes(x = group, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size = 17)) + geom_point() + 
        annotate("segment", x=min(wss$group), xend=trunc(h.life), y=wss[trunc(h.life),'wss'], yend=wss[trunc(h.life),'wss'], colour = "red", linetype = "longdash") +
        annotate("text", x=1, y=wss[trunc(h.life),'wss'],
                 label=paste('h.life'), vjust=-.5, size=6, colour='red') +
        annotate("text", x=trunc(h.life), y=wss[trunc(h.life),'wss'],
                 label=paste(trunc(h.life)), vjust=-.5, size=6, colour='red')
      
      if (!file.exists("./plots/")) dir.create("./plots/")
        ggsave(paste0("./plots/Kmeans_clusterAnalysis_",projName, ".png"), dpi = 300)
      dev.off()
    } else {
      ggplot(wss, aes(x = group, y = wss)) + geom_line() + labs(x = "Number of Clusters", y = "Within groups sum of squares") + theme(text = element_text(size = 17)) + geom_point() +
        annotate("rect", xmin=which(wss$RatioChange==min(wss$RatioChange[-1]))-1, 
                 xmax=which(wss$RatioChange==min(wss$RatioChange[-1])), 
                 ymin=min(wss$wss), ymax=wss[which(wss$RatioChange==min(wss$RatioChange[-1])),'wss'], alpha=.3, fill = "red") +
        annotate("text", x=which(wss$RatioChange==min(wss$RatioChange[-1]))-1, 
                 y=wss[which(wss$RatioChange==min(wss$RatioChange[-1])),'wss'], 
                 label=paste(which(wss$RatioChange==min(wss$RatioChange[-1]))-1), vjust=-1.5, size=5, colour='red')
      if (!file.exists("./plots/")) dir.create("./plots/")
      ggsave(paste0("./plots/Kmeans_clusterAnalysis_",projName, ".png"), dpi = 300)
      dev.off()
      }
    
  }
  
  if (randomforest) {
    cat("Processing randomForest \n")
    #Generating random points:
    if (is.null(random.pt)) {
      random.pt <- as.integer(0.01 * ncell(envLayer))
    }
    cat("Generating", random.pt,"random points \n")
    set.seed(seed)
    sp.pts <- spsample(studyArea, random.pt, type = "random")
    
    #Extracting raster values for each point
    sp.pts$vals <- extract(envLayer, sp.pts)
    
    #removing NA
    if (any(is.na(sp.pts@data))) {
      sp.pts <- sp.pts[-which(is.na(sp.pts@data)),]
    }
    
    #creating object with the raster value to run kmeans
    mydata <- as.data.frame(sp.pts@data)
    
    # Not used for this prupouse
    mydata <- scale(mydata) # standardize variables
    
    # K-Means Cluster Analysis
    fit <- kmeans(
      mydata, ngroup, algorithm = c("Hartigan-Wong",
                                    "Lloyd",
                                    "Forgy",
                                    "MacQueen")[1], nstart = 10)
    
    # Saving statistical results
    if(save.fit) save(fit, file=paste0(folder, projName, '_RF_','.RData'))
    
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
    
    if (save.shp) {
      cat("converting and saving to Polygons \n")
      # Saving Vector output (Polygon)
      rf_segmentation <-
        rasterToPolygons(
          rf_segmentation, fun = NULL, n = 4, na.rm = TRUE, digits = 12, dissolve = FALSE)
      writeOGR(
        rf_segmentation, dsn = paste0(folder), layer = paste0('rf_segmentation_',projName), driver =
          "ESRI Shapefile", overwrite_layer = TRUE
      )
    }
    if (save.raster){
      cat("Saving raster result \n")
      #Saving RASTER output
      if (!file.exists("./rasters/")) dir.create("./rasters/")
      writeRaster(
        rf_segmentation, filename = paste0(folder, './rasters/rf_segmentation_',projName,'.tif'), overwrite = TRUE
      )
    }
    cat("randomForest segmentation done. \n")
  }
  if (Kmeans) {
    cat("Processing K-means \n")
    mydata <- getValues(envLayer)
    
    # Identificando pixel nao NA
    a <- which(!is.na(mydata))
    
    # Creating data to Kmeans analysis
    mydata <- na.omit(values(envLayer))
    
    # Not used for this prupouse
    mydata <- scale(mydata) # standardize variables
    
    # K-Means Cluster Analysis
    set.seed(seed)
    fit <-
      kmeans(
        mydata, ngroup, iter.max = 50, algorithm = c("Hartigan-Wong",
                                      "Lloyd",
                                      "Forgy",
                                      "MacQueen")[1], nstart=10)
    
    # Saving statistical results
    if(save.fit) save(fit, file=paste0(folder, projName, '_Kmeans_', '.RData'))
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
          km_SegmentationRaster, dsn = paste0(folder), layer = paste0('km_segmentation_',projName), driver = "ESRI Shapefile", overwrite_layer = TRUE)
      } 
    if (save.raster) {
        #Saving RASTER output
        cat("Saving raster result \n")
      if (!file.exists("./rasters/")) dir.create("./rasters/")
        writeRaster(
          km_SegmentationRaster, filename = paste0(folder, './rasters/km_segmentation_',projName,'.tif'), overwrite = TRUE)
        }
    cat("Kmeans segmentation done. \n")
  }
  if (fuzzy.cluster){
    cat("Processing C-means (fuzzy K-means) \n")
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
                    rate.par = 0.3, weights = 1, control = list())
    
    # Saving statistical results
    if(save.fit) save(fuzzy, file=paste0(folder, projName, '_Fuzzy_','.RData'))
    
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
        fuzzy_SegmentationRaster, dsn = paste0(folder), layer = paste0('fuzzy_segmentation_',projName), driver = "ESRI Shapefile", overwrite_layer = TRUE)
    } 
    if (save.raster) {
      #Saving RASTER output
      cat("Saving raster result \n")
      if (!file.exists("./rasters/")) dir.create("./rasters/")
      writeRaster(
        fuzzy_SegmentationRaster, filename = paste0(folder, './rasters/fuzzy_segmentation_',projName,'.tif'), overwrite = TRUE)
    }
    cat("Fuzzy segmentation done. \n")
  }
  result <- list()
  if (exists('rf_segmentation')) result<-c(result, rf_segmentation)
  if (exists('km_SegmentationRaster')) result <-c(result, km_SegmentationRaster)
  if (exists('fuzzy_SegmentationRaster')) result<-c(result, fuzzy_SegmentationRaster)
  return(result)
  }