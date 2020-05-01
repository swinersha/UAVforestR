
#' Takes an image and extracts a spatialPolygonsDataframe of the tree crown locations.
#'
#' @param imagery A raster image (usually pre blurred)
#' @param im_sobel A raster image the same size as imagery but containing its sobel convolution
#' @param THRESHSeed The crown minimum height as a proportion (which is relative to
#' maximum crown height)
#' @param THRESHCrown The crown minimum height as a proportion (which is relative to
#' mean crown height, which is adjusted on the fly during processing)
#' @param tau the allometric percentile used to limit tree crown growth
#' @param spectT The minimum height to be considered in the segmentation
#' @param SOBELstr The sobel stength used to detect edges; this is a little cryptic
#' to interpret. It is set relative to the height of the tree and essentially
#' encodes the edge strength needed to define an edge. Larger values of SOBELstr
#' will make the edge detect less sensitive. Values between 0 and 5 are usually a good
#' starting point.
#' itcIMG_fast
#' @param lm.searchwin  The search window for finding maxima. If this is set to NULL the window will scale
#' allometrically according to tree height. Otherwise it can be set to an odd integer which will
#' be fixed for the detection of all maxima, irrespective of tree height.
#' @param pypath The path to gdal_polygonize.py on the system; if this is set to NULL, the function
#' raster::rasterToPolygons will be used but please not that this is very slow for large images.
#' @return A spatialPolygonsDataframe of tree crown locations
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-08-17

##########################################################
# A modification of ITCsegment package by Michele Dalponte

# Changes made by Tom Swinfield for the RSPB under funding from the Cambridge Conservation Initiative
# Collaborative Fund
##########################################################

# This is a reverse engineering of ITC segment topdown to first find all the local maxima and
# then to work through these from tallest to shortest using the fast algorithm.

# Extracts just the locations of the tree maxima:
itcIMG_fast <- function (imagery,
                         im_sobel,
                         THRESHSeed,
                         THRESHCrown,
                         lut,
                         tau,
                         specT,
                         SOBELstr,
                         lm.searchwin = NULL,
                         pypath
                         )
{
  # Extracts the image data as a matrix:
  img <- matrix(dim(imagery)[2], dim(imagery)[1], data = imagery[,], byrow = FALSE)
  Sobel <- matrix(dim(imagery)[2], dim(imagery)[1], data = im_sobel[,], byrow = FALSE)
  img <- img[1:dim(imagery)[2], dim(imagery)[1]:1]
  Sobel <- Sobel[1:dim(imagery)[2], dim(imagery)[1]:1]
  img[is.na(img)] <- 0 # Sets any nas to 0s.
  img[img < specT] <- 0 # Any values beneath the minimum height are set to 0s.

  # Finds the local maxima:
  coordSeeds <-
    detect.maxima(img, res(imagery)[1], lm.searchwin = lm.searchwin, lut = lut, tau = tau)
  coordCrown <- segment.crowns(
    x = img,
    x.Sobel = Sobel,
    coordSeeds,
    THRESHSeed = THRESHSeed,
    THRESHCrown = THRESHCrown,
    SOBELstr = SOBELstr,
    scale = res(imagery)[1],
    lut = lut,
    tau = tau
  )

  crowns_sp<-crowns_to_spatial(coordCrown, imagery, pypath, convex.hull = FALSE)
  maxima_sp<-maxima_to_spatial(coordCrown, imagery, pypath)
  crowns_sp<-crown_max_combine(crowns_sp, maxima_sp)

  maxima_sp<<-maxima_sp

  return(crowns_sp)
}
