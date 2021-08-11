

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


