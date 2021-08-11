#'
#' Fills the working terrain file using kNN, iterating
#'   using the specified known WSE.  Elevations outside 
#'   +/- 1 ft of this elevation are considred valid 
#'   bathymetry/terrain grid cells and will be used
#'   in the interpolation.
#'


### Potential future work:
# could add some criterion to the iteration algorithm, but I don't 
#   think that iterating more on the kNN will be a problem unless it
#   is exceedingly time-consuming.  Otherwise, something like this
#   might be necessarry:

# Testing if predictions are converging or not from previous prediction
#  Criterion is that 75% of points are 
# if(sum(abs(rXYZ[rXYZ$fillPts, "z"] - (wse-0.5))/sum(rXYZ$fillPts)) > 0.75)

library(rgdal)
library(raster)
library(FNN)
library(sp)
library(tidyverse)

### Config -------------------------------------------------------------------

# relative to data/working directory
workingDir <- "R06"

# Elevations outside of +/- 1 ft. of this elevation will be used as valid
#   interpolation points
wse <- 392.

# Bathymetry raster file name in workingDir
rFileName <- "Terrain6.tif"

# Polygon of grid cells to fill overlain on terrainFileName,
#   file in workingDir
fillShpFileName <- "R06_fill.shp"

# Max. number of kNN iterations to attempt before stopping algorithm
maxkNNIts <- 500

# k-nearest neighbor value
k = 200


# Options are: localAverage, kNN
regMethod <- "localAverage"

### Objects ------------------------------------------------------------------

# A new type of matrix, where cell queries outside of the range of the matrix 
#   (i.e, row/cell < 1 or greater than nrow/ncol) would return a default
#   value.  This would be ideal for approximation algorithms so that 
#   you don't need to worry about boundary condition issues.
#
# Currently, a matrix request outside of the cell boundary retuns this:
# > mat[-1:1,-1:1]
# Error in mat[-1:1, -1:1] : 
#   only 0's may be mixed with negative subscripts
#
# Ideally, this would return either a matrix or vector of points, where
#   the invalid points are filled with the default value.
#
# Use the setClass method:
# https://stat.ethz.ch/R-manual/R-devel/library/methods/html/setClass.html

### Functions ----------------------------------------------------------------

#' For each point indicated that needs to be smoothed, a simple approximation
#'   method will be used to approximate the point with the neraby points
#'
smoothMatrix <- function(m, smoothPts){
  
}


### Load ---------------------------------------------------------------------

# Bathy raster
rawR <- r <- 
  sprintf("data/working/%s/%s", workingDir, rFileName) %>%
  raster()

# Fill polygon
rawFill <- 
  sprintf("data/working/%s/%s", workingDir, fillShpFileName) %>%
  readOGR()

### Main ---------------------------------------------------------------------


# TESTING 
# 

# Crop the terrain to just a subset to speed up calc testing
testExtent <- extent(-5720000, -5700000, 9190000, 9200000) # xmin, xmax, ymin, ymax
rawR <- r <- crop(x = rawR, y = testExtent)
plot(r)

#
# END TESTING



# Extract initial XYZ, coords from raw raster [X,Y,Z]
rXYZ <- rawR %>% coordinates() %>% as.data.frame() %>% 
  mutate(z = values(r),
         i = rep(1:ncol(r), times = nrow(r)),
         j = rep(1:nrow(r), each = ncol(r)))
names(rXYZ) <- c("x","y","z", "i", "j")

# Extracting polygon XY [X,Y]
pXY <- rawFill@polygons[[1]]@Polygons[[1]]@coords %>% as.data.frame()
names(pXY) <- c("x","y")

# Selecting those points of raster within fill polygon and indicating
#   them as points to iteratively compute kNN
rXYZ$fillPts <- 
  sp::point.in.polygon(point.x = rXYZ$x, point.y = rXYZ$y,
                       pol.x = pXY$x, pol.y = pXY$y) %>%
  as.logical()

# Adding another layer of logic to the fillPts layer, restricting them to
#   within +/- 1 ft of the 'wse' defined in config
rXYZ$fillPts <- rXYZ$z < wse
rXYZ$fillPts[is.na(rXYZ$fillPts)] <- F

# Iteratively filling in previous k-NN prediction 
#   of nearest neighbors for points
for(i in 1:maxkNNIts){
  cat(sprintf("\n\t%03g\t%s", i, Sys.time()))
  
  
  if(regMethod == "localAverage"){
    
    m <- as.matrix(rawR)
    rawR <- smoothMatrix(m = m, smoothPts = rXYZ$fillPts)
    
    
  }
  
  
  if(regMethod == "kNN"){
    # On first iteration, just use known valid points
    if(i == 1){
      rXYZ$z[rXYZ$fillPts] <- 
        knn.reg(train =  rXYZ[!rXYZ$fillPts, c("x","y")],
                test = rXYZ[rXYZ$fillPts, c("x","y")],
                y = rXYZ[!rXYZ$fillPts, "z"],
                k = k)$pred
    }else{
      # After furst iteration, use all points
      rXYZ$z[rXYZ$fillPts] <- 
        knn.reg(train =  rXYZ[rXYZ$fillPts, c("x","y")],
                test = rXYZ[rXYZ$fillPts, c("x","y")],
                y = rXYZ[rXYZ$fillPts, "z"],
                k = k)$pred
    }
  }
}

# Reassigning data to raster
setValues(r, rXYZ$z)
plot(r)

# Save raster to new location
saveFileName <- sprintf("data/product/%s/%s", workingDir, rFileName)
if(!dir.exists(dirname(saveFileName))) dir.create(dirname(saveFileName))
raster::writeRaster(x = r, filename = saveFileName, overwrite=T)
