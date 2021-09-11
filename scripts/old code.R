

### Old Code - Finding best WSE to subtract souding depths from


#' computeSoundingError <- function(w, bathyXYZ, soundingXYZ, k=3){
#'   #'
#'   #' Function to compute the differences between the bathymetric
#'   #'   points and sounding points in 3D Eulidian distances
#'   #'   via kNN.  The differences are returned as a sum
#'   #'
#'   library(FNN)
#'   soundingXYZ$z <- w - soundingXYZ$z
#'   get.knnx(data = bathyXYZ ,query = soundingXYZ, k = k)$nn.dist %>% 
#'     rowSums() %>% sum()
#' }
#' 
#' # Minimizing both the difference between:
#' #   1) the sounding and bathymetry points
#' # by gradually adjusting the wse in brute force
#' bestW <- 0
#' minE <- 1E100
#' for(w in seq(wse-3,wse+3, by=0.1)){
#'   # Apply new WSE, retrieving error
#'   err <- computeSoundingError(w, bathyXYZ, soundingXYZ, k=100)
#'   if(err < minE) {
#'     bestW <- w
#'     minE <- err
#'   }
#' }


### Old Code - formatting ----------------------------------------------------


# 
# # Determining which points are water surface points by seeing if
# #   any of the differences in elevation of surrounding points is less than
# #   a theshold value
# # Initial kNN, excluding first nearest neighbor, which is current point
# rawKNN <- FNN::knnx.index(data = rawXYZ, query = rawXYZ, k = 9)
# # Matrix of values of nearest neighbors
# rawValMat <- matrix(data = rawXYZ$z[rawKNN[,2:9]],
#                     nrow = nrow(rawKNN), ncol = 8)
# elevMat <- matrix(data = rep(rawXYZ$z[rawKNN[,2:9]],each = 8),
#                   nrow = nrow(rawKNN), ncol = 8, byrow = T)
# diffMat <- abs(elevMat - rawValMat)
# maxDiffs <- apply(X = diffMat, MARGIN = 1, FUN = max)
# 
# wseIndex <- interpIndex <- maxDiffs < 0.001
# 
# bathyIndex <- rawXYZ$z < wse-wseBuffer & !wseIndex
# 
# bathyXYZ <- rawXYZ[bathyIndex,] %>% as.data.frame()
# # Establishing index points below assumed flat WSE
# # wseIndex <- interpIndex <- rawXYZ[,3] <= wse+wseBuffer & 
# # rawXYZ[,3] >= wse-wseBuffer
# wseXYZ <- rawXYZ[wseIndex,] %>% as.data.frame()
# soundingXYZ <- rawSounding %>% as.data.frame() %>%
#   dplyr::select("coords.x1","coords.x2","coords.x3")
# names(soundingXYZ) <- names(wseXYZ) <- names(bathyXYZ) <- 
#   names(modXYZ) <- names(rawXYZ) <- c("x","y","z")
# 

