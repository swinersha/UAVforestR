#-----------------------------------------------
# 
# A function that finds the polygon with the highest overlap for a given polygon. 
#
# It was written specifically to match polygons of automatically segmented tree crowns
# with polygons of manually segmented tree crowns.
#
# Tom Swinfield
# 17-02-16
#
#-----------------------------------------------

#auto_trees=ut; manual_trees=manual; buffer_by=60; verbose='on'
#i<-1
#auto_trees<-ut; manual_trees<-manual[1:5,]; buffer_by=60

crown_overlap<-function(auto_trees, manual_trees, buffer_by, verbose='off'){
  out<-matrix(0, nrow=nrow(manual_trees), ncol=4)
  auto_trees$id<-1:nrow(auto_trees)
  sum(sapply(auto_trees@polygons, function(x) x@area)<0.001)
  for(i in 1:nrow(manual_trees)){
#    print(i)
    poly<-manual_trees[i,] # selects the current polygon
    poly_area<-poly@polygons[[1]]@Polygons[[1]]@area # the area of the polygon
    poly_buffer<-gBuffer(poly, width=buffer_by) # makes a buffer around it
    cropped_trees<-raster::crop(auto_trees, poly_buffer) # crops the auto trees to the buffer
    cropped_trees <- gBuffer(cropped_trees, byid=TRUE, width=0) # deals with self-intersection 
    for(j in 1:nrow(cropped_trees)){
#      print(j)
      overlap<-gIntersection(poly, cropped_trees[j,]) # extracts the intersecting area.
      if(!is.null(overlap)){
        auto_area<-cropped_trees@polygons[[j]]@Polygons[[1]]@area
        overlap_area<-overlap@polygons[[1]]@Polygons[[1]]@area # the area of the intersection
        overlap_percent<-overlap_area/poly_area # the percentage area
        size_ratio<-auto_area/poly_area # the percentage area
        if(verbose=='on')
          cat('j: ', j, 'overlap: ', overlap_percent, '\n')
        if(overlap_percent>out[i,3]){ # stores the result
          out[i,1]<-i
          out[i,2]<-cropped_trees$id[j]
          out[i,3]<-overlap_percent
          out[i,4]<-size_ratio
        }
      }
    }
    if(verbose=='on')
      cat('out: ', out[i,], '\n')
  }
  out[out[,2]==0,2]<-NA
  # Loads these as additional columns for the manual trees:
  manual_trees$id_auto_tree<-out[,2]
  manual_trees$overlap_auto_tree<-out[,3]
  manual_trees$size_ratio<-out[,4]
  return(manual_trees)
}    