#'
#' Establishing global variables
#'

library(raster)
library(yaml)
library(purrr)
library(tidyverse)
library(kriging)

### Config -------------------------------------------------------------------

# Config location
cnfg <- read_yaml("config.yaml")

# Save directory relative to project directory
saveDir <- "interp_raster"

### Functions ----------------------------------------------------------------

loadTerrain <- function(fileName){
  #' Loads raster dataset using global 'cnfg' object and relative path 
  #'   to the file name
  library(raster)
  raster(x = file.path(getwd(), cnfg$dataDir, fileName))
}


getCLine <- function(r=NA, cLineFileName=NA){
  #'
  #' Compute the centerline of a river based upon the approximation of the
  #'   thalweg
  #'
  
  return(NA)
}

interpTerr <- function(rawTerr, thalwegCline){
  #'
  #' Interpolate terrain using kriging, radiating from the centerline
  #'   outward to the floodpains.
  #'
  return(NA)
}

saveRasters <- function(saveDir, saveTer){
  #'
  #' Save to TIFF raster
  #'
  return(NA)
}

splitDfEvenly <- function(inDf, nVals=1000000){
  #' Even split a dataframe into a list of dataframes, each
  #'   with nVals rows
  split(inDf, ceiling(seq_along(1:nrow(inDf))/nVals))
}

krigPts <- function(inDf){
  #' Krigs a matrix of points where the 
  #'   second dimension is x, y, and z
  kriging(inDf[,1], inDf[, 2], inDf[,3], model="gaussian")
}

renameRawPts <- function(inCols){
  inCols[grepl("Terrain",inCols)] <- "z"
  inCols
}

### Main ---------------------------------------------------------------------

terFileNames <- map_chr(.x = cnfg$data,.f = function(x)x$fileName)
bFileNames <- tools::file_path_sans_ext(basename(terFileNames))

# Load terrain datasets
rawTer <- terFileNames %>% map(loadTerrain)

# Convert to points
rawPts <- rawTer %>% map(rasterToPoints)

rawDf <- rawPts %>% map(as.data.frame) %>% 
  map(rename_with, renameRawPts) %>%
  map2_df(bFileNames, ~ mutate(.x, ID = .y)) 

# Create new column indicating where to split each dataframe

# Set thalweg centerline
# thalwegCLine <- getCLine(ter)

# Interpolate along thalweg for NA grid cells
# saveTer <- interpTerr(rawTer, thalwegCLine)


# Krig and replace values where z is less than specified elev

kPts <- rawPts %>% map(krigPts)

# Write to save directory
saveRasters(saveDir, saveTer)

