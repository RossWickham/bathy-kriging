#'
#' Class definition for a new type of flexible matrix 'flexMat' that
#'   allows you toquery for index ranges outside of the matrix, returning
#'   a default value for those items that don't exist.  The default
#'   value is defined when creating the flexible matrix object.
#'


library(pryr)


flexMat <- structure(matrix(), class = "flexMat")
f <- function(x) UseMethod("f")

#'
#' Defining the flexMat class to be initialized from an existing matrix
#'
flexMat <- function(m, defaultValue=NA){
  
  
  ## Get the environment for this
  ## instance of the function.
  thisEnv <- environment()
  
  m <- m
  
  ## Create the list used to represent an
  ## object for this class
  out <- list(
    # m = m,
    # defaultValue = defaultValue,
    # ## Define the environment where this list is defined so
    # ## that I can refer to it later.
    # thisEnv = thisEnv,
    m = get("m",thisEnv),
    
    ## Returns a matrix or vector of cell values
    #   corresponding to the desired cell range
    getCells = function(xmin, xmax, ymin, ymax, asVector=TRUE,
                        defaultValue = get("defaultValue", thisEnv))
    {
      vals <- matrix(0, nrow = xmax - xmin + 1, ncol = ymax - ymin + 1)
      
      # Fill columns on left side of matrix
      if(xmin < 1) vals[1:-xmin,] <- defaultValue
      # Fill rows on top of matrix
      if(ymin < 1) vals[,1:-ymin] <- defaultValue
      # Fill columns on right side of matrix
      if(xmax > nrow(vals)) vals[,(ymax-nrow(vals)):nrow(vals)] <- defaultValue
      # Fill rows on bottom of matrix
      if(ymax < ncol(vals)) vals[(xmax-ncol(vals)):ncol(vals),] <- defaultValue
      
      if(asVector) vals <- as.numeric(vals)
      vals
    }
    
  )
  
  ## Define the value of the list within the current environment.
  assign("m",out,envir=thisEnv)
  class(out) <- append(class(out), "flexMat")
  out
}

# Testing
fM <- flexMat(matrix(1:100, nrow=10))

class(fM)
str(fM)

fM$getCells(-2,4,-1,4, asVector = F)
