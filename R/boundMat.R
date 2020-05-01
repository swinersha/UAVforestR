#' Creates a matrix for the space around local maxima to store information about
#' the reasons for boundary detection.
#'
#' @description  The matrix is scaled by tree height in order
#' to maintain computational speed.
#' @param r the row position within the image (used as a spatial offset term)
#' @param k the column position within the image (used as a spatial offset term)
#' @param scale the one dimensional pixel size / resolution of the image
#' @param tau the allometric percentile used to limit tree crown growth
#' @return A list containing 1. an empty boundMat matrix; 2. the row offset term;
#' 3. the column offset term; 4 the maximum pixel distance to the crown boundary (
#' this is the same as the distance from the centre of the matrix to its edge).
#' @export
#' @author Tom Swinfield
#' @details Created 17-08-17

boundMat_set <- function(r, k, z, scale, lut, tau) {
  pix2edge <- edge.dist(z, scale = scale, lut = lut, tau = tau)
  focalWinSize <- pix2edge * 2 + 3
  boundMat <- matrix(0, focalWinSize, focalWinSize)
  boundMat <- list(boundMat, r, k, pix2edge)
  return(boundMat)
}

#' Update boundMat
#'
#' @description A table used to update the boundMat matrix according to the values in the focal
#' pixels being parsed by itcIMG_fast. It is basically a helper function.
#' @param boundMat A boundMat object (created by boundMat_set)
#' @param filData The boundary conditions of the pixels being parsed by itcIMG_fast
#' @param fil.dists The distance to each pixel being parsed from the tree maximum / seed.
#' @param i An integer encoding the tree identity
#' @param THRESHSeed The crown minimum height as a proportion (which is relative to
#' maximum crown height)
#' @param THRESHCrown The crown minimum height as a proportion (which is relative to
#' mean crown height, which is adjusted on the fly during processing)
#' @param rvSobel The sobel strength, as scaled by the tree height, needed to set a boundary position.
#' @param crownrad.THRESH The maximum crown radius set by the allometric relationship with height.
#' @return An updated boundMat object
#' @export
#' @author Tom Swinfield
#' @details Created 17-08-17

boundMat_edit <-
  function(boundMat,
           filData,
           fil.dists,
           i,
           rvSeed,
           rvSobel,
           THRESHCrown,
           THRESHSeed,
           crownrad.THRESH) {
  coords <- filData[, 1:2]
  coords[, 1] <- coords[, 1] - boundMat[[2]] + boundMat[[4]] + 1
  coords[, 2] <- coords[, 2] - boundMat[[3]] + boundMat[[4]] + 1

  coords_esc<<-coords

  # This is in reverse order so that the last statements (which are the most important) dominate.
  boundMat[[1]][coords] <- 7

  boundMat[[1]][coords[filData[, 3] > (rvSeed + (rvSeed * 0.05)), , drop=FALSE]] <-
    6 # pixel taller than tree maximum
  boundMat[[1]][coords[filData[, 5] >= rvSobel, , drop=FALSE]] <-
    5 # Pixel greater than edge threshold
  boundMat[[1]][coords[filData[, 3] <= (rvSeed * THRESHCrown), , drop=FALSE]] <-
    4 # pixel less than mean height threshold
  boundMat[[1]][coords[filData[, 3] <= (rvSeed * THRESHSeed), , drop=FALSE]] <-
    3 # pixel less than max height threshold
  boundMat[[1]][coords[!(filData[, 4] %in% c(0, i)), , drop=FALSE]] <-
    2 # pixel in another tree crown
  boundMat[[1]][coords[fil.dists >= crownrad.THRESH, , drop=FALSE]] <-
    1 # pixel greater than allometric threshold
  return(boundMat)
}

#' Calculate the total number of crown boundary pixels defined according
#' to each boundary condition.
#'
#' @param boundMat A boundMat object (created by boundMat_set)
#' @return An integer vector the same length as the number of boundary conditions. At present this
#' is hardcoded to 6. The integer 7 is used to code pixels in the crown but not part of the
#' boundary.
#' @export
#' @author Tom Swinfield
#' @details Created 17-08-17
boundMat_summarise<-function(boundMat){
  bm_sum<-sapply(1:6, function(i) sum(boundMat[[1]] == i))
  return(bm_sum)
}

#' View boundMat matrix focussed on its centre.
#'
#' @description A helper function to view the part of boundMat where the action is happening
#' @param boundMat A boundMat object (created by boundMat_set)
#' @param size The size of the matrix to print.
#' @return The cropped matrix focussed on the crown maximum / seed.
#' @export
#' @author Tom Swinfield
#' @details Created 17-08-17
#'
boundMat_view<-function(boundMat, size = 9){
  if((size %% 2) == 0)
    size<-size + 1
  index<-(boundMat[[4]] - (size - 1)/2 ) : (boundMat[[4]] + (size - 1)/2 + 1)

  bm_view<-boundMat[[1]][index, index]
  return(bm_view)
}

#' Calculate the average sobel edge strength at the crown boundary.
#'
#' @description A helper function to calculate the average of the sobel edge image for the crown
#' boundary. This function only averages pixels which were segmented based upon either
#' the height or sobel boundary conditions.
#'
#' @param boundMat A boundMat object (created by boundMat_set)
#' @param img_sobel The sobel image
#' @return The average sobel edge strength as a number.
#' @export
#' @author Tom Swinfield
#' @details Created 17-08-17
#'
boundMat_sobel_mean<-function(boundMat, img_sobel){
  coords<-lapply(3:5, function(i) which(boundMat[[1]] == i, arr.ind = TRUE))
  coords<-do.call(rbind, coords)
  sobel_mean<-mean(img_sobel[coords], na.rm = TRUE)
  return(sobel_mean)
}



