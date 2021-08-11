#'
#' Example of config creation for example dataset
#' 


library(yaml)
out <- 
  list(data = list(
    list(name="R06",
         fileName="R06/Terrain6.tif"),
    list(name="R07",
         fileName="R07/Terrain7.tif"),
    list(name="R08",
         fileName="R08/Terrain8.tif")
  ),
  dataDir = "_data"
  )
yaml::write_yaml(out,"config.yaml")
