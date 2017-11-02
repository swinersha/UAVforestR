
#' Converts the coordCrown object created by itcIMG_fast into a spatialPolygonsDataframe
#' of segmented tree crowns
#'
#' @param coordCrown A list of pixel positions attributed to each tree by id, produced by
#' itcIMG_fast
#' @param x The original image
#' @param pypath The path to gdal_polygonize.py on the system.
#' itcIMG_fast
#' @param convex.hull Logical, indicating whether or not the spatial polygons should just be
#' encircled with a convex hull or not.
#' @return A spatialPolygonsDataframe of tree crown locations
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-08-17

crowns_to_spatial <- function(coordCrown, x, pypath, convex.hull=FALSE) {
  Crowns_img <- x

  Crowns <- matrix(ncol(x), nrow(x), data = 0, byrow = FALSE)
  Crowns[coordCrown[[1]][, 1:2]] <- coordCrown[[1]][, 3]
  Crowns_img[] <-
    as.numeric(Crowns[1:dim(Crowns)[1], dim(Crowns)[2]:1],
               byrow = TRUE)
  Crowns_img[Crowns_img == 0] <- NA # Excludes non maxima pixels.

  Crowns_img <- methods::as(Crowns_img, "SpatialGridDataFrame")
  Crowns_img <- raster::raster(Crowns_img, layer = 1)

  if (!is.null(pypath))
    Crowns_sp <-
    gdal_polygonizeR(
      Crowns_img,
      outshape = NULL,
      gdalformat = 'ESRI Shapefile',
      pypath = pypath,
      readpoly = TRUE,
      quiet = TRUE
    )
  else
    Crowns_sp <- raster::rasterToPolygons(Crowns_img, fun = , dissolve = TRUE)

  # Extracts the mean and max tree heights:
  names(Crowns_sp@data) <- "tree"
  Crowns_sp <- Crowns_sp[Crowns_sp$tree != 0, ]
  Crowns_sp@data <-
    cbind(Crowns_sp@data, coordCrown[[2]][Crowns_sp$tree, ])

  Crowns_sp$A <-
    unlist(lapply(Crowns_sp@polygons, function(x)
      methods::slot(x, "area")))
  Crowns_sp$R <- sqrt(Crowns_sp$A / pi)

  # Create convex hulls around the polygons:
  if (convex.hull) {
    Crowns_sp <-
      rgeos::gBuffer(Crowns_sp,
                     width = -res(imagery)[1] / 2,
                     byid = T)
    ITCcv <- rgeos::gConvexHull(Crowns_sp, byid = T)
    Crowns_sp <-
      sp::SpatialPolygonsDataFrame(ITCcv, data = Crowns_sp@data, match.ID = F)
  }
  return(Crowns_sp)
}


#' Converts the coordCrown object created by itcIMG_fast into a spatialPolygonsDataframe
#' of tree maxima locations
#'
#' @param coordCrown A list of pixel positions attributed to each tree by id, produced by
#' itcIMG_fast
#' @param x The original image
#' itcIMG_fast
#' @param pypath The path to gdal_polygonize.py on the system.
#' @return A spatialPolygonsDataframe of tree crown locations
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-08-17
maxima_to_spatial <- function(coordCrown, x, pypath) {
  Maxima_img <- x

  Maxima <- matrix(ncol(x), nrow(x), data = 0, byrow = FALSE)
  Maxima[as.matrix(coordCrown[[2]][, 1:2])] <- 1:nrow(coordCrown[[2]])
  Maxima_img[] <-
    as.numeric(Maxima[1:dim(Maxima)[1], dim(Maxima)[2]:1],
               byrow = TRUE)
  Maxima_img[Maxima_img == 0] <- NA # Excludes non maxima pixels.

  Maxima_img <- methods::as(Maxima_img, "SpatialGridDataFrame")
  Maxima_img <- raster::raster(Maxima_img, layer = 1)
  if (!is.null(pypath))
    Maxima_sp <-
    gdal_polygonizeR(
      Maxima_img,
      outshape = NULL,
      gdalformat = 'ESRI Shapefile',
      pypath = pypath,
      readpoly = TRUE,
      quiet = TRUE
    )
  else
    Maxima_sp <-
    raster::rasterToPolygons(Maxima_img, fun = , dissolve = TRUE)
  return(Maxima_sp)
}


#' Combines the crown maxima locations with the crown polygons.
#'
#' @description The crown polygons spatialPolygonsDataframe is returned
#' but with new columns max.x and max.y describing the spatial coordinates of the maxima.
#'
#' @param crowns The crowns spatialPolygonsDataframe
#' @param maxima The maxima spatialPolygonsDataframe
#' @return A spatialPolygonsDataframe of tree crown polygons
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-08-17
crown_max_combine<-function(crowns, maxima) {
  match.id <- match(crowns$tree, maxima$DN)
  maxima_coords <- coordinates(maxima)[match.id,]
  colnames(maxima_coords) <- c("max.x", "max.y")
  crowns@data<-cbind(crowns@data, maxima_coords)
  return(crowns)
}
