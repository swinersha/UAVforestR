
#' Run gdal_polygonize from R
#'
#' @description A fast function for creating a SpatialPolygonsDataFrame from a raster:
# https://johnbaumgartner.wordpress.com/2012/07/26/getting-rasters-into-shape-from-r/
#' @param x A raster image
#' @param outshape The name of the shapefile if saving
#' @param gdal_format The format of the shapefile to use; by default ESRI shapefile format
#' @param pypath The path to gdal_polygonize.py on the system
#' @param readpoly Logical, encoding whether or not the output should be returned or not. If
#' not NULL is returned.
#' @param quiet Logical indicating whether or not a verbose output should be creaeted.
#' @return A spatialPolygonsDataframe
#' @export
#' @author John Baumgartner
#' @details
#'
#' Created 17-08-17


gdal_polygonizeR <- function(x, outshape=NULL, gdalformat = 'ESRI Shapefile',
                             pypath=NULL, readpoly=TRUE, quiet=TRUE) {
  if (isTRUE(readpoly)) require(rgdal)
  if (is.null(pypath)) {
    pypath <- Sys.which('gdal_polygonize.py')
  }
  #if (!file.exists(pypath)) stop("Can't find gdal_polygonize.py on your system.")
  owd <- getwd()
  on.exit(setwd(owd)) # very handy little function
  setwd(dirname(pypath))
  if (!is.null(outshape)) {
    outshape <- sub('\\.shp$', '', outshape)
    f.exists <- file.exists(paste(outshape, c('shp', 'shx', 'dbf'), sep='.'))
    if (any(f.exists))
      stop(sprintf('File already exists: %s',
                   toString(paste(outshape, c('shp', 'shx', 'dbf'),
                                  sep='.')[f.exists])), call.=FALSE)
  } else outshape <- tempfile()
  if (is(x, 'Raster')) {
    require(raster)
    writeRaster(x, {f <- tempfile(fileext='.tif')}) # Also handy bit of code
    rastpath <- normalizePath(f)
  } else if (is.character(x)) {
    rastpath <- normalizePath(x)
  } else stop('x must be a file path (character string), or a Raster object.')
  #system2('python', args=(sprintf('"%1$s" "%2$s" -f "%3$s" "%4$s.shp"',
  #                                pypath, rastpath, gdalformat, outshape)))
  system2(command = pypath, args=(sprintf('"%1$s" -f "%2$s" "%3$s.shp"',
                                  rastpath, gdalformat, outshape)))

  if (isTRUE(readpoly)) {
    shp <- readOGR(dirname(outshape), layer = basename(outshape), verbose=!quiet)
    return(shp)
  }
  return(NULL)
}
