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

getShpLineCoords <- function(shp){
  coordinates(shp)[[1]][[1]] %>%
    as.data.frame() %>%
    setNames(c("x","y"))
}

### Testing ------------------------------------------------------------------


# TESTING 
# 
if(F){
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

crsObj <- crs(rawRaster)

xyz <- initializeXYZ(rawRaster, rawSounding, wse)

# Ensure centerline same projection as raster
rawCL <- spTransform(rawCL, crsObj)

rasterFromXYZ(xyz[, c("x","y","z")]) %>%
  plot()
rawCL %>% lines()

### Bathymetric Lateral Interpolation ----------------------------------------


#
# TODO
# - apply smoothing on interpolated cells after finish compute
# - add more vertices in centerline
#


if(T){
  
  #'
  #' Create the 'X' matrix of points from xyz data.frame
  #'
  getX <- function(xyz){
    matrix(data = c(xyz$x, xyz$y), nrow = 2, byrow = 2)
  }
  
  #'
  #' Create a rotational matrix given a vertical slope
  #'   i.e., the horizontal:vertical ration is 1:s
  #'
  getR <- function(s){
    theta = atan(s)
    matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow = 2)
  }
  
  #'
  #' Create the shift matrix to be applied to the 'X' matrix
  #'
  getS <- function(xShift, yShift, nPts){
    matrix(c(rep(-xShift, nPts), rep(-yShift, nPts)), nrow = 2, byrow = T)
  }
  
  #'
  #' Compute current centerline slope from index value and centerline
  #'
  computeCLSlope <- function(ix, xlCL){
    with(xlCL,
         1/2*(
           (y[ix+1] - y[ix])/(x[ix+1] - x[ix]) + 
             (y[ix] - y[ix-1])/(x[ix] - x[ix-1])
         )
    )
  }
  
  #'
  #' Create a matrix structured where top row is x-values,
  #'   second row is y-values.
  #'
  getXYMatrix <- function(xyz){
    matrix(c(xyz$x, xyz$y), nrow = 2, byrow = T) 
  }
  
  #'
  #' Find a set of equally spaced points on a line given two vertices
  #'   defining the end points 
  #'
  #' @param xy data.frame with two columns for x and y, and two
  #'             rows defining
  #' @param d  numeric, atomic defining spacing of points to
  #'             to define along line
  #' @return   data.frame with two columns for x and y, 
  #'             defining an sequence of evenly-spaced (d) points
  #'             between the extents of the input xy data.frame
  #'
  getEqualSpacePtsOnLine <- function(xy, d){
    # Compute a linear model fit from xy
    lmFit <- 
      xy %>%
      lm(y ~ x, .)
    # Get the slope of the line in radians
    theta <- 
      lmFit %>%
      coefficients() %>%
      last() %>%
      atan()
    # Compute the horizontal (x) distance given the 
    #   hypotenuse (d).  cos(theta) = d/x_dist =>
    x_dist <- d*cos(theta) %>% abs()
    # Generate a sequence from start to end of x-values in xy
    if(diff(range(xy$x)) < x_dist){
      return(
        # Return mean of x and y values if equal spacing in x-direction
        #   is less than the extent defined in xy input data.frame
        xy %>% sapply(mean) %>% t() %>% as.data.frame()
      )
    }else{
      # Get a sequence between extents of xy's x values, starting at x_dist, 
      #   and compute the y-values using the linear model
      xValues <- seq(from = x_dist, to = max(xy$x), by = x_dist)
      yValues <- predict.lm(object = lmFit, newdata = data.frame(x = xValues))
      return(
        data.frame(
          x = xValues,
          y = yValues
        )
      )
    }
  }
  
  #'
  #' Generate a sequence of points at which to apply the 
  #'   inpterolation algorithm.  Returns a dataframe
  #'   with columns
  #'
  #' @param clXY data.frame of the centerline from which to apply
  #'              interpolation with columns for x and y values
  #' @param a  (Optional) numeric defining the rotation angle to increment
  #'             when interpolating at corners
  #' @return a new data.frame with columns defining interpolation points with
  #'           columns for x-values, y-values, slope, and theta (radians) 
  getInterpolationPOints(clXY, a = 0.1){
    
    # Compute data.frame with slopes and theta values associated 
    #   with each line segment for use later
    
    # Between each vertex, generate a sequence of points
    
    # - Extract two-row xy data.frame, 
    # - pull the points
    # - compute a slope
    # - compute a theta = atan(slope)
    #
    
    # At each corner (non-last or non-first vertice), create a set
    #   of points rotating 'a' radians from the slope of the previous
    #   line segment to the slope of the next.
    
    # - use slope and theta values from above computed segment data.frame
    # - determine the angle of rotation from the smallest angle between
    #     previous and next slopes defining vertice
    # - generate a sequence of angles between the two
    # - assign angles to new data.frame using fixed x-y values of corner
    
    NULL
  }
  
  xlCL <- getShpLineCoords(rawCL)
  
  l_buff <- 500
  t_buff <- 4000
  
  for(i in 2:(nrow(xlCL)-1)){
    
    
    # Reserve a dummy set of points for manipulation
    mXYZ <- xyz
    
    # Compute centerline slope
    s <- computeCLSlope(i, xlCL)
    
    # Get the rotational matrix
    R <- getR(s = s)
    
    # Get the shift matrix
    S <- getS(xShift = xlCL$x[i], yShift = xlCL$y[i], nPts = nrow(mXYZ))
    
    # Shift and rotate
    mXYZ[, c("x","y")] <-t(R %*% (getXYMatrix(mXYZ) + S))
    
    # Select points within l_buff and t_buff rectangle about newly defined 
    #   origin for interpolation
    ptsInInterpRectanlge <- 
      mXYZ$x > -l_buff & 
      mXYZ$x < l_buff &
      mXYZ$y > -t_buff &
      mXYZ$y < t_buff
    if(sum(ptsInInterpRectanlge) == 0) next
    
    sXYZ <- mXYZ[ptsInInterpRectanlge,]
    
    # plot(sXYZ[, c("y","z")])
    
    # Add checks for 1) some non-NA points, 2) no NA points
    if(sum(is.na(sXYZ$z)) < 50) next
    if(!any(is.na(sXYZ$z))) next
    
    # Use good points in this subset of points to define the bad points
    knnQuery <- 
      get.knnx(data = sXYZ[!is.na(sXYZ$z), c("x","y")],
               query = sXYZ[is.na(sXYZ$z), c("x","y")],
               k = 3)
    
    zValues <- sXYZ$z[!is.na(sXYZ$z)]
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
              1/x^2
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
    sXYZ$z[is.na(sXYZ$z)] <- 
      ((kValues * kWeights) *
         kNormalize) %>%
      apply(MARGIN = 1,
            FUN = sum)
    
    # Assign back to master
    xyz$z[ptsInInterpRectanlge] <- sXYZ$z
    
    
  }
  
  rasterFromXYZ(xyz[, c("x","y","z")]) %>%
    plot()
  
  cat(sprintf("\n\nWriting raster to here:\n\t%s", saveFileName))
  xyz[, c("x","y","z")] %>%
    rasterFromXYZ(crs=crsObj) %>%
    writeRaster(filename = saveFileName, overwrite=T)
}
