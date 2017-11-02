
#' Caclulate the the performance of itfIMGfast parameters against a training set of manually segemented crowns.
#'
#' @description Finds the automatically segmented tree crown with the highest overlap with a given manually
#' segemented polygon tree crown. The function was written specifically to match polygons of automatically segmented tree crowns
#' with polygons of manually segmented tree crowns.
#' @param auto_trees A set of tree crowns automatically produced by itcIMG_fast
#' @param manual_trees A spatialPolygonsDataframe of manually segmented tree crowns
#' @param buffer_by The area to buffer around the focal crown when looking for matching crowns
#' @param verbose Logical indicating whether or not the progress should be printed
#' @return A matrix of coordinates of tree maxima locations in pixels and tree heights
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-02-16

# auto_trees=ut; manual_trees=manual; buffer_by=60; verbose='on'
# i<-51; j<-1
# auto_trees<-ut; manual_trees<-manual[1:5,]; buffer_by=60

crown_overlap<-function(auto_trees, manual_trees, buffer_by, verbose='off'){
  # out<-matrix(0, nrow=nrow(manual_trees), ncol=4)
  out<-matrix(0, nrow=nrow(manual_trees), ncol=6)
  out[,1]<-1:nrow(manual_trees)
  out[,3]<-1

  auto_trees$id<-1:nrow(auto_trees)
  sum(sapply(auto_trees@polygons, function(x) x@area)<0.001)
  for(i in 1:nrow(manual_trees)){
    i_esc<<-i
#    print(i)
    poly<-manual_trees[i,] # selects the current polygon
    poly_area<-poly@polygons[[1]]@Polygons[[1]]@area # the area of the polygon
    poly_buffer<-gBuffer(poly, width=buffer_by) # makes a buffer around it
    #cropped_trees<-raster::crop(auto_trees, poly_buffer) # crops the auto trees to the buffer
    cropped_trees<-raster::crop(auto_trees, poly) # crops the auto trees to any that intersect with poly
    cropped_trees<-auto_trees[auto_trees$id %in% cropped_trees$id,]
    cropped_trees <- gBuffer(cropped_trees, byid=TRUE, width=0) # deals with self-intersection
    if(!is.null(cropped_trees)){
      for (j in 1:nrow(cropped_trees)) {
        #      print(j)
        overlap <-
          gIntersection(poly, cropped_trees[j, ]) # extracts the intersecting area.
        if (!is.null(overlap)) {
          auto_area <- cropped_trees@polygons[[j]]@Polygons[[1]]@area
          overlap_area <-
            overlap@polygons[[1]]@Polygons[[1]]@area # the area of the intersection

          #overseg <- (auto_area - overlap_area) / poly_area
          overseg <- 1-(overlap_area / auto_area) # false positive rate
          underseg <- (poly_area - overlap_area) / poly_area # false negative rate
          # overlap_percent<-overlap_area/poly_area # the percentage area
          size_ratio<-auto_area/poly_area # the percentage area

          if (out[i, 3] == 1) {
            out[i, 2] <- cropped_trees$id[j]
            out[i, 3] <- 0
            out[i, 4] <- overseg
            out[i, 5] <- underseg
            out[i, 6] <- size_ratio
          }
          else if ((overseg + underseg) < sum(out[i, 4:5])) {
            # stores the result
            out[i, 2] <- cropped_trees$id[j]
            out[i, 4] <- overseg
            out[i, 5] <- underseg
            out[i, 6] <- size_ratio
          }
          # if(verbose=='on')
          #   cat('j: ', j, 'overlap: ', overlap_percent, '\n')
          # if(overlap_percent>out[i,3]){ # stores the result
          #   out[i,1]<-i
          #   out[i,2]<-cropped_trees$id[j]
          #   out[i,3]<-overlap_percent
          #   out[i,4]<-size_ratio
          # }
        }
      }
      if (verbose == 'on')
        cat('out: ', out[i, ], '\n')
    }
  }
  out[out[,2]==0,2]<-NA
  # Loads these as additional columns for the manual trees:
  manual_trees$id_auto_tree<-out[,2]
  # manual_trees$overlap_auto_tree <- out[, 3]
  # manual_trees$size_ratio <- out[, 4]
  manual_trees$tree_match <- out[, 3]
  manual_trees$overseg <- out[, 4]
  manual_trees$underseg <- out[, 5]
  manual_trees$cost <- rowSums(out[,4:5])/2 + out[,3]
  manual_trees$size_ratio<-out[,6]

  if (verbose == 'on')
    cat('Cost: ', manual_trees$cost, '\n')

  return(manual_trees)
}
