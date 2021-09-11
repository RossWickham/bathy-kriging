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
wse <- 392
# wse <- 440

# Buffer distance at which to include points for interpolation
wseBuffer <- 2

# Sounding polypoint shapefile
shpFileName <- "data/raw/Harbour_Sounding_point.shp"

# k-nearest neighobors to use for estimation
k = 50

# Number of smoothing k-NN iterations
nSmoothingIts <- 4

### Load ----------------------------------------------------------------------

rawRaster <- r <- rFileName %>% raster()

# Load shapefile, projecto to common crs
rawSounding <- s <- shpFileName %>% readOGR() %>% spTransform(CRSobj = crs(r))

# Load centerline
rawCL <- readOGR(dsn = "data\\working\\R06\\R06_cl.shp")

### Functions ----------------------------------------------------------------

#'
#' Adding new columns for row and column corresponding to a matrix from a 
#'   gridded dataframe with columns for 'x' and 'y'spatial coordinates 
#' @input  inDf dataframe with named 'x' and 'y' columns defining a uniform grid
#' @output dataframe with new columns for 'row' and 'col' indicating the
#'         numeric 1:n values for rows and columns.
addRowColumnToGriddedDf <- function(inDf){
  rows <- sort(unique(inDf$x))
  columns <- sort(unique(inDf$y))
  inDf$row <- match(inDf$x, rows)
  inDf$col <- match(inDf$y, columns)
  inDf
}


#'
#' Defines new columns in a dataframe for points defining the bathymetric
#'   sounding locations and the water surface elevation.  A buffer elevation
#'   is used.
#' 
#'
computeBathyAndWSEIndices <-
  function(xyz, wse=NULL, wseBuffer = 1, bathyXYZ = NULL){
    # Initial kNN eval for all points
    rawKNN <- 
      FNN::knnx.index(data = xyz[,c("x","y")],
                      query = xyz[,c("x","y")],
                      k = 9)
    # Matrix of values of nearest neighbors
    rawValMat <- matrix(data = xyz$z[rawKNN[,2:9]],
                        nrow = nrow(rawKNN), ncol = 8)
    # Elevation matrix of current point
    elevMat <- matrix(data = rep(xyz$z[rawKNN[,2:9]],each = 8),
                      nrow = nrow(rawKNN), ncol = 8, byrow = T)
    # Computing difference between each point and surrounding points
    diffMat <- abs(elevMat - rawValMat)
    # Computing max diference between each point and surrounding points in grid
    maxDiffs <- apply(X = diffMat, MARGIN = 1, FUN = max)
    
    # Determining the WSE index based on given threshold
    xyz$wseIndex <- xyz$interpIndex <- maxDiffs < 0.001
    # If defined, use WSE and WSE buffer to define additiona WSE index 
    #   locations
    if(!is.null(wse))
      xyz$wseIndex <- 
      xyz$interpIndex <- 
      xyz$wseIndex | 
      (xyz$z >= wse - wseBuffer &  xyz$z <= wse + wseBuffer)
    # Determining the points corresponding to the bathymetry data
    xyz$bathyIndex <- xyz$z < wse-wseBuffer & !xyz$wseIndex
    xyz
  }

#'
#' Determines the compute order by iteratively going through
#'   each 
#' @input xyz             dataframe with gridded 'row' and 'column' 
#'          column names
#' @input interpIndexColName column corresponding to the 
defineInterpComputeOrder <- function(xyz, interpIndexColName){
  nComputeCells <- sum(xyz$interpIndex)
  computeCellOrder <- 1 # initialize compute cell order
  xyz$computeOrder <- NA
  xyz$dummyInterpIndex <- xyz[,interpIndexColName]
  pass <- 1
  
  
  ### TODO
  # May want to improve the interpolation method by using a cell adjacent to 
  #   previous "good" cell as a starting point for the next interpolation point
  # Also, why does it iteration on passes without finding a vali point?  
  #   Seems like it should find a good point eventually
  
  
  # Iterating through all row-column combination, searching for cells 
  #   that have NA values and filling in a dummy column to track if a cell
  #   has already been tallied as a compute cell.
  while(sum(xyz$dummyInterpIndex) > 0){
    # Randomly walk through all potential cells to use for interpolation
    for(n in sample(which(xyz$dummyInterpIndex), replace = F)){
      i <- xyz$row[n]
      j <- xyz$col[n]
      surroundingValues <- 
        getSurroundingCellValues(xyz, rowIndex = i, colIndex = j)
      # xyz$dummyInterpIndex[xyz$xyz >= i-1 & xyz$row <= i+1 & 
      #            xyz$col >= j-1 & xyz$col <= j+1]
      cat(sprintf("\npass = %5g\tn = %6g", pass, n))
      if(any(!is.na(surroundingValues))){
        cat(sprintf("\ti = %g\tj = %g, # remaining = %g",
                    i, j, nComputeCells - computeCellOrder))
        xyz$dummyInterpIndex[n] <- F
        xyz$computeOrder[n] <- computeCellOrder
        computeCellOrder <- computeCellOrder + 1
      }
    }
    pass <- pass + 1
  }
  cat("\nDone defining cell compute order")
  xyz$dummyInterpIndex <- NULL
  xyz
}


#'
#' Returns the indices of the surrounding cells
#'
isSurroundingCell <- 
  function(xyz, computeOrder=NULL, rowIndex=NULL, colIndex = NULL, ix=NULL){
    if(!is.null(computeOrder))
      rowValues <- xyz[match(x = 1, table = xyz$computeOrder),]
    if(!is.null(rowIndex) & !is.null(colIndex))
      rowValues <- xyz[xyz$row == rowIndex & xyz$col == colIndex,]
    if(!is.null(ix)) rowValues <- xyz[ix,]
    xyz$row >= rowValues$row-1 & 
      xyz$row <= rowValues$row+1 &
      xyz$col >= rowValues$col-1 &
      xyz$col <= rowValues$col+1
  }

#'
#' Pull values of surrounding cells given either:
#'   1) the compute order or
#'   2) the row and column index values
#'
getSurroundingCellValues <- 
  function(xyz, computeOrder=NULL, rowIndex=NULL, colIndex = NULL, ix=NULL){
    indices <- isSurroundingCell(xyz, computeOrder, rowIndex, colIndex, ix)
    xyz$z[indices]
  }

#'
#' Iteratively loop through compute order, taking average
#'   of all nearby non-NA cells until a threshold difference in
#'   error is acheived between current and previous cell values. 
#'
interpolateCells <- function(xyz, maxErr = 0.1, maxIts = 1000){
  # Reducing to just cells with a compute order, preserving index
  #   of cells needing to be 
  
  # Initialize
  xyz$prevZ <- xyz$newZ <- xyz$z
  it <- 1
  currentMaxError <- 1E10
  maxComputeOrder <- max(xyz$computeOrder, na.rm = T)
  # Order dataframe by compute order, row number is now compute order
  xyz <- xyz[order(xyz$computeOrder, na.last = T),]
  while(currentMaxError >= maxErr){
    cat(sprintf("\nIteration:\t%6g\tError:\t%.1f", it, currentMaxError))
    for(n in 1:nrow(xyz)){
      if(is.na(xyz$computeOrder[n])) break
      if(n %% 1000 == 0)
        cat(sprintf("\n\tCompute order:\t%10g/%10g", n, maxComputeOrder))
      surroundingValues <- getSurroundingCellValues(xyz, n)
      xyz$newZ[n] <- mean(na.omit(surroundingValues))
    }
    maxErrEstimate <- max(abs(xyz$prevZ - xyz$newZ), na.rm=T)
    if(maxErrEstimate < currentMaxError & it != 1)
      currentMaxError <- maxErrEstimate
    cat(sprintf("\nDone with iteration %6g.  Maximum absolute error:\t%.1f",
                it, maxErrEstimate))
    it <- it + 1
    if(it > maxIts | currentMaxError < maxErr) break
    xyz$prevZ <- xyz$newZ # recylce latest prediction
  }
  xyz$prevZ <- NULL
  xyz
}

#'
#' Adds a new column to a dataframe indicating the
#'   number NA values in surrounding cells
#'
computeNumberSurroundingNAandMinValue <- function(xyz){
  xyz$naCount <- 0
  xyz$minSurroundingValue <- NA
  for(n in which(xyz$interpIndex)){
    surroundingCellValues <- 
      getSurroundingCellValues(xyz,
                               rowIndex = xyz$row[n],
                               colIndex = xyz$col[n])
    surroundingCellValues <- 
      surroundingCellValues[!is.infinite(surroundingCellValues) &
                              !is.na(surroundingCellValues)]
    xyz$naCount[n] <- sum(is.na(surroundingCellValues))
    xyz$minSurroundingValue[n] <- 
      min(surroundingCellValues, na.rm=T)
  }
  xyz
}

knnSmoothXYZ <- function(xyz, k=10){
  xyz$z[!is.na(xyz$z)] <- 
    knn.reg(train = xyz[!is.na(xyz$z) ,c("x","y")],
            test  = xyz[!is.na(xyz$z),c("x","y")],
            y = xyz$z,
            k=k)$pred
  xyz
}

initializeXYZ <- 
  function(rawRaster, rawSounding, wse){
    
    # Extracting xyz location of bathymetric points
    bathyXYZ <-
      rawSounding %>% 
      as.data.frame() %>% 
      dplyr::select(c("coords.x1","coords.x2", "coords.x3")) %>%
      rename(x = coords.x1, y = coords.x2, z = coords.x3)
    
    # Characterize points using 'wse', extracting xyz, renaming columns
    xyz <-
      rawRaster %>%
      rasterToPoints() %>%
      as.data.frame() %>% 
      rename(z=Terrain6) %>%
      addRowColumnToGriddedDf()  %>%
      mutate(
        wseIndex = z >= wse - wseBuffer & z <= wse + wseBuffer,
        interpIndex = wseIndex,
        bathyIndex = F
      )
    
    # Finding the nearest raster points' index corresponding
    #   to the bathymetric points, update those points with z-values
    #   and excluding them from interpolation
    knnQuery <- 
      knnx.index(data = xyz[c("x","y")],
                 query = bathyXYZ[c("x","y")],
                 k = 1) %>% 
      as.numeric()
    
    # Use tidyverse syntax here to create/manipulte column info
    xyz$z[knnQuery] <- xyz$z[knnQuery] - bathyXYZ$z
    xyz$interpIndex[knnQuery] <- xyz$wseIndex[knnQuery] <- F
    xyz$bathyIndex[knnQuery] <- T
    
    # Setting interpolation cells to NA
    xyz$z[xyz$interpIndex] <- NA
    
    xyz[,c("x","y","z")] %>%
      rasterFromXYZ() %>%
      plot()
    points(bathyXYZ[, c("x","y")])
    
    xyz
  }

### Testing ------------------------------------------------------------------


# TESTING 
# 
if(T){
  # Crop the terrain to just a subset to speed up calc testing
  testExtent <- extent(-5720000, -5700000, 9190000, 9192000) # xmin, xmax, ymin, ymax
  rawRaster <- r <- crop(x = rawRaster, y = testExtent)
  
  
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
  
  plot(r)
  points(s)
}
#
# END TESTING

# Main ------------------------------------------------------------------------

# set.seed(100)

crsString <- crs(rawRaster)

xyz <- initializeXYZ(rawRaster, rawSounding, wse)

rasterFromXYZ(xyz[, c("x","y","z")]) %>%
  plot()


### Bathymetric Lateral Interpolation ----------------------------------------
if(T){



}
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