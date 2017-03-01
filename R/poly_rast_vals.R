#------------------------------------------------
# A simple function to extract raster values given a spatial polygon
# object. The extracted values are then summarised given a function.

# Tom Swinfield
# 17-02-16

#------------------------------------------------

# A function enabling raster value extraction and summarisation: 
poly_rast_vals<-function(poly, img, FUN){
  # extracts the data for each polygon from the raster:
  poly_vals = raster::extract(img,poly)
  # applies a function to these:
  fun_vals<-sapply(poly_vals, FUN=FUN)
  return(fun_vals)
}