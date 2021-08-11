#'
#' Code to fill in a partial bathymetry dataset by combining:
#'   1. The original bathy rater
#'   2. A raster of valid points in the raster
#'   3. A polypoint of depth soundings
#'   
#' The general process for computing the raster is as follows:
#'   - poly
#'
#'



# spatial
library(rgeos)
library(raster)
library(sp)
library(rgdal)
library(gstat)

#
library(rgl)
library(purrr)
library(tidyverse)
library(FNN)

### Config -------------------------------------------------------------------

# Bathy raster
rFileName <- "data/raw/R06/Terrain6.tif"

# Save file
# saveFileName <- "data/product/R06/Terrain6_knn.tif"
saveFileName <- "D:\\projects\\h_h\\snake_temperature\\data\\raster\\Terrain6_knn.tif"

# The water surface elevation to be used for determining which points in the
#   raster are valid bathymetry points, and which need to be interpolated
# wse <- 392
wse <- 440

# Buffer distance at which to include points for interpolation
wseBuffer <- 1

# Sounding polypoint shapefile
shpFileName <- "data/raw/Harbour_Sounding_point.shp"

# k-nearest neighobors to use for estimation
k = 5

# Number of smoothing k-NN iterations
nSmoothingIts <- 4

### Load ----------------------------------------------------------------------

rawR <- r <- rFileName %>% raster()

# Load shapefile, projecto to common crs
rawSounding <- s <- shpFileName %>% readOGR() %>% spTransform(CRSobj = crs(r))


### Testing ------------------------------------------------------------------


# TESTING 
# 
if(T){
  # Crop the terrain to just a subset to speed up calc testing
  testExtent <- extent(-5720000, -5700000, 9190000, 9200000) # xmin, xmax, ymin, ymax
  rawR <- r <- crop(x = rawR, y = testExtent)
  
  
  gClip <- function(shp, bb){
    library(sp)
    if("matrix" %in% class(bb)){
      ex <- extent(as.vector(t(bb)))
    }else{
      ex <- extent(bb)
    }
    crds_poly <- ex %>% coordinates()
    crds_shp <- coordinates(shp)
    inBbox <-
      point.in.polygon(point.x = crds_shp[,1], point.y = crds_shp[,2],
                       pol.x = crds_poly[,1], pol.y = crds_poly[,2]) %>%
      as.logical()
    shp[inBbox,]
  }
  
  rawSounding <- s <- gClip(rawSounding, bbox(r))
  
  # Saving shapefile in proper crs to get *.prj file
  rgdal::writeOGR(obj = rawSounding,
                  driver = "ESRI Shapefile",
                  layer = basename(tools::file_path_sans_ext(saveFileName)),
                  dsn = paste0(tools::file_path_sans_ext(saveFileName),".shp"))
  
  plot(r)
  points(s)
}
#
# END TESTING




# Main ------------------------------------------------------------------------

# Characterize points using 'wse', extracting xyz, renaming columns
rawXYZ <- modXYZ <- rawR %>% rasterToPoints() %>% as.data.frame()
names(modXYZ) <- names(rawXYZ) <- c("x","y","z")

# Determining which points are water surface points by seeing if
#   any of the differences in elevation of surrounding points is less than
#   a theshold value
# Initial kNN, excluding first nearest neighbor, which is current point
rawKNN <- FNN::knnx.index(data = rawXYZ, query = rawXYZ, k = 9)
# Matrix of values of nearest neighbors
rawValMat <- matrix(data = rawXYZ$z[rawKNN[,2:9]],
                    nrow = nrow(rawKNN), ncol = 8)
elevMat <- matrix(data = rep(rawXYZ$z[rawKNN[,2:9]],each = 8),
                  nrow = nrow(rawKNN), ncol = 8, byrow = T)
diffMat <- abs(elevMat - rawValMat)
maxDiffs <- apply(X = diffMat, MARGIN = 1, FUN = max)

wseIndex <- interpIndex <- maxDiffs < 0.001

bathyIndex <- rawXYZ$z < wse-wseBuffer & !wseIndex

bathyXYZ <- rawXYZ[bathyIndex,] %>% as.data.frame()
# Establishing index points below assumed flat WSE
# wseIndex <- interpIndex <- rawXYZ[,3] <= wse+wseBuffer & 
  # rawXYZ[,3] >= wse-wseBuffer
wseXYZ <- rawXYZ[wseIndex,] %>% as.data.frame()
soundingXYZ <- rawSounding %>% as.data.frame() %>%
  dplyr::select("coords.x1","coords.x2","coords.x3")
names(soundingXYZ) <- names(wseXYZ) <- names(bathyXYZ) <- 
  names(modXYZ) <- names(rawXYZ) <- c("x","y","z")


### Slow Averaging -----------------------------------------------------------
# Slowly iterating through the cells selected for interpolation that have
#   at least one valued (i.e., non-NA) cell adjacent to them, while also 
#   keeping track of the order of approximation.  Iteratively pass though
#   the order of cells to approximat until a desired difference in change 
#   is achieved using a while loop.


# Add row and column identifiers to raw matrix





### Kriging ------------------------------------------------------------------
# Try kriging, following example here:
# https://rpubs.com/nabilabd/118172
coordinates(modXYZ) = ~x+y
lzn.vgm <- variogram(object = log(z)~x+y, data = modXYZ[!interpIndex,]) # calculates sample variogram values 
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Mat", 900, 1)) # fit model
plot(lzn.vgm, lzn.fit)

# Kriging
lzn.kriged <- krige(log(z)~x+y, bathyXYZ, wseXYZ, model=lzn.fit)

modXYZ$z[wseIndex] <- exp(as.numeric(lzn.kriged@data$var1.pred))

# If this looks good, then might want to smooth using knn
rgl::plot3d(rawXYZ[!bathyIndex & !wseIndex,])
points3d(modXYZ, col="red")


### kNN ----------------------------------------------------------------------


# Using the best WSE to fit the data
# Finding the nearest raster point's index in bathyXYZ
#   corresponding to each soundingXYZ
knnQuery <- knnx.index(data = rawXYZ[c("x","y")],
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


rgl::plot3d(rawXYZ[!bathyIndex & !wseIndex,])
points3d(modXYZ, col="red")
points3d(bathyXYZ, col="green", size=2)

rasterFromXYZ(modXYZ, crs=crs(r)) %>%
  writeRaster(filename = saveFileName, overwrite=T)
