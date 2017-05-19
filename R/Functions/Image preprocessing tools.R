#------------------------------------------------
# A couple of wrappers for running the sobel edge detection and blurring as part
# of the pre-processing for ITC segment

# Tom Swinfield
# 17-02-16

#------------------------------------------------


sobel_edge<-function(imagery){
  # Calculating edge strength using a Sobel edge detector:
  kernX<-matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3)
  kernY<-matrix(c(-1,-2,-1,0,0,0,1,2,1), 3, 3, byrow=TRUE)
  im_sobel <- raster::focal(imagery, w = matrix(1, 3, 3), 
        fun = function(x) {
            y<-sqrt(sum(x* kernX)^2 + sum(x* kernY)^2)
            if(is.na(y))
              y<-0
            return(y)
            })
  return(im_sobel)
}

# Blurs an image with a standard gaussian blur:
blur<-function(imagery){
  imagery <- raster::focal(imagery, w = matrix(1, 3, 3), 
          fun = function(x) {
            mean(x, na.rm = T)
         })
  return(imagery)
}
