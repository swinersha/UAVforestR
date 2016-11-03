# imgPoly_match takes a vector of file paths to suitable images and checks their spatial extents 
# against a polygon to see if they are included, at least in part within the polygon.

#img.path=chm.tiles; poly=comps
#x<-img.path[1]
imgPoly_match<-function(img.path, poly)
{
  # works through the file paths to iamges one by one to check for the best match:  
  # this could be made way more efficient with a bit of thought.
  tile.test<-lapply(img.path, function(x)
  {
    img<-raster::brick(x)
    if(sp::proj4string(poly)!=sp::proj4string(img))
      poly<-sp::spTransform(poly, crs(img))
    ei1 <- as(raster::extent(img), "SpatialPolygons")
    ei2 <- as(raster::extent(poly), "SpatialPolygons")
    rgeos::gContainsProperly(ei1, ei2) # test whether the raster contains the polygon or not.
  }
  )
  img.path[unlist(tile.test)]
}

#--------------------------------------------------
# 161102
# This is an updated version that is more efficient and powerful in the way it can be specified
# Takes a directory, a regular expression pattern, and a spatial polygons dataframe
# it returns a list of strings which are file paths to the images that the polygons 
# are contained within:

#dir<-chm.path; pattern=paste('Chunk', chunk, 'CHM', '.*.tif$')
#dir<-om.path; pattern=paste('Chunk', chunk, 'OM', '.*.tif$')
#list.files(dir)


imgPoly_match2<-function(dir, pattern=NULL, poly)
{
  
  if(is.null(pattern))
     files<-file.path(dir, list.files(dir))
  else
     files<-file.path(dir, list.files(dir, pattern=pattern))
     
     rasters<-lapply(files, raster) # Loads the rasters
     
     rasters_crs<-sapply(rasters, function(X) raster::crs(X)) # Find the CRS of all the rasters
     rasters_extent<-sapply(rasters, raster::extent) # Find the extent of all the rasters
     u.crs<-unique(rasters_crs)
     # Check that we are not dealing with multiple coordinate reference systems:
     if(length(u.crs)>1)
       stop("Error: multiple CRSs in use.")
     poly<-sp::spTransform(poly, u.crs[[1]]) # Tranform the polygons to match the images
     
     # works through the file paths to iamges one by one to check for the best match:  
     # this could be made way more efficient with a bit of thought.
     poly_tiles<-lapply(1:nrow(poly), function(i)
     {
       poly_sub<-poly[i,]
       ei1 <- as(raster::extent(poly_sub), "SpatialPolygons") # find the extent of the polygon
       y<-rasters[[1]]
       tile.test<-sapply(rasters_extent, function(y)
       {
         ei2 <- as(y, "SpatialPolygons") # the extent of the image
         rgeos::gContainsProperly(ei2, ei1) # test whether the raster contains the polygon or not.
       }
       )
       files[unlist(tile.test)] # index from the list of filepaths to the images
     }
     )
     return(poly_tiles)
}

