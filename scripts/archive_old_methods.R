#'
#' Archive of previous attempts to interpolate
#'
#'
#'


### Bottom Up ----------------------------------------------------------------
# Compute from the lowest points up, selecting cells with NA 
#   surrounding values


# Smooth non-NA cells first
# xyz <- 
#   xyz %>%
#   knnSmoothXYZ()


## TODO
# Add more inteligent interpolation
#  only use cells within x distance for interpolation
#  inverse distance method
#  interpolate all cells at once via FNN, assign NA if no close matches 
#    within x specified distance for interpolation

if(F){
  
  maxDist <- 100  # feet
  
  it <- 1
  maxIts <- 1000
  while(T){
    
    interpIndices <- xyz$interpIndex
    
    kDist <- knnx.dist(data =  xyz[!interpIndices, c("x","y")],
                       query = xyz[, c("x","y")],
                       k=k)[,-1]
    minDists <-
      apply(X = kDist,
            MARGIN = 1,
            FUN = min)
    
    # Determine compute order
    interpIndices <- 
      is.na(xyz$z) |
      xyz$z == 0 &
      minDists < maxDist
    
    if(sum(interpIndices) == 0){
      cat("\nNo more indexes to interpolate, escaping")
      break
    }
    
    # Get kNN indices and distances for cell locations to compute 
    #   in this iteration
    zValues <- xyz$z[!xyz$interpIndex]
    knnQuery <- 
      get.knnx(data = xyz[!interpIndices, c("x","y")],
               query = xyz[interpIndices, c("x","y")],
               k=k)
    # Pull the values of nearest neighbors for each index
    kValues <- 
      apply(X = knnQuery$nn.index,
            MARGIN = 1,
            FUN = function(x){
              zValues[x]
            }) %>%
      t()
    
    # Convert the distances to weights
    kWeights <- 
      apply(knnQuery$nn.dist,
            MARGIN = 1,
            FUN = function(x){
              out <- 1/x^2
              out[x > maxDist] <- 0
              out
              # 1/x^2
            }) %>%
      t()
    # Normalizing matrix to use as numerator
    kNormalize <- 
      apply(kWeights,
            MARGIN = 1, 
            FUN = function(x){
              1/rep(sum(x),length(x))
            }) %>%
      t()
    
    
    # Computing estimates
    xyz$z[interpIndices] <- 
      ((kValues * kWeights) *
         kNormalize) %>%
      apply(MARGIN = 1,
            FUN = sum)
    
    
    interpIndices[apply(kWeights,1,max) ==0] <- F
    
    # Update indices that need to be interpolated
    xyz$interpIndex[interpIndices] <- F
    
    if(sum(xyz$interpIndex)==0) break
    cat(sprintf("\nits = %5g\t# cells remaining:\t%10g", it, sum(xyz$interpIndex)))
    it <- it + 1
    if(it > maxIts) break
  }
  
  
  rasterFromXYZ(xyz[, c("x","y","z")]) %>%
    plot()
  
  if(file.exists(saveFileName)) file.remove(saveFileName)
  
  xyz[, c("x","y","z")] %>%
    rasterFromXYZ(crs=crsString) %>%
    writeRaster(filename = saveFileName, overwrite=T)
}

### Slow Averaging -----------------------------------------------------------
# Slowly iterating through the cells selected for interpolation that have
#   at least one valued (i.e., non-NA) cell adjacent to them, while also 
#   keeping track of the order of approximation.  Iteratively pass though
#   the order of cells to approximat until a desired difference in change 
#   is achieved using a while loop.

if(F){
  
  maxErr <- 0.1
  maxIt <- 1000
  xyz$err <- 0
  # Make an initial pass through the raster locations to determine
  #   the order in which to interpolate the gridded points
  
  
  # Iterating through each point in the WSE index
  while(abs(max(xyz$err)) > maxErr & it < maxIt){
    for(it in 1:maxIt){
      for(k  in which(xyz$wseIndex)){
        
      }
    }
  }
  
}


### Kriging ------------------------------------------------------------------
# Try kriging, following example here:
# https://rpubs.com/nabilabd/118172

if(F){
  coordinates(modXYZ) = ~x+y
  # calculates sample variogram values
  lzn.vgm <- variogram(object = log(z)~x+y, data = modXYZ[!interpIndex,])  
  lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Mat", 900, 1)) # fit model
  plot(lzn.vgm, lzn.fit)
  
  # Kriging
  lzn.kriged <- krige(log(z)~x+y, bathyXYZ, wseXYZ, model=lzn.fit)
  
  modXYZ$z[wseIndex] <- exp(as.numeric(lzn.kriged@data$var1.pred))
  
  # If this looks good, then might want to smooth using knn
  rgl::plot3d(xyz[!bathyIndex & !wseIndex,])
  points3d(modXYZ, col="red")
  
}
### kNN ----------------------------------------------------------------------

if(F){
  
  # Using the best WSE to fit the data
  # Finding the nearest raster point's index in bathyXYZ
  #   corresponding to each soundingXYZ
  knnQuery <- knnx.index(data = xyz[c("x","y")],
                         query = soundingXYZ[c("x","y")],
                         k=1) %>% as.numeric()
  
  # Updating which points need interpolation, excluding sounding points
  interpIndex[knnQuery] <- F 
  
  
  ## TODO 
  # interpolate a WSE raster based on the sounding points that
  #   overlap with bathy XSecs.  WSE = bathy + sounding depth
  # Interpolate between known WSEs to fill in the raster
  
  # Assigning those raster values to the corresponding
  #   sounding depth
  # modXYZ$z[knnQuery] <- bestW - soundingXYZ$z
  
  # Subtracting soundings from WSE to estimate bathymetric elevations
  modXYZ$z[knnQuery] <- modXYZ$z[knnQuery] - soundingXYZ$z
  
  # Iteratively interpolating, starting by first only
  #   informing interpolation with valid terrain/bathymetry points
  modXYZ$z[interpIndex] <- 
    knn.reg(train = modXYZ[!interpIndex, c("x","y")],
            test = modXYZ[interpIndex, c("x","y")],
            y = modXYZ[!interpIndex, "z"] %>% as.data.frame(),
            k = k)$pred
  
  # Smoothing bathymetry points
  for(i in 1:nSmoothingIts){
    cat(sprintf("\n%g", i))
    modXYZ$z[interpIndex] <-
      knn.reg(train = modXYZ[!interpIndex, c("x","y")],
              test = modXYZ[interpIndex, c("x","y")],
              y = modXYZ[!interpIndex, "z"] %>% as.data.frame(),
              k = k)$pred
  }
  
  
  rgl::plot3d(xyz[!bathyIndex & !wseIndex,])
  points3d(modXYZ, col="red")
  points3d(bathyXYZ, col="green", size=2)
  
  rasterFromXYZ(modXYZ, crs=crs(r)) %>%
    writeRaster(filename = saveFileName, overwrite=T)
}