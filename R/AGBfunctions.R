

cover<-function(x, h_thresh, na.rm){
  na.ind<-is.na(x)
  if(any(na.ind))
    return(NA)
  x<-x[!na.ind]
  n<-length(x)
  y<-(sum(x>=h_thresh))/n
  return(y)
}

cover_n<-function(x, h_thresh, na.rm){
  na.ind<-is.na(x)
  na.n<-sum(na.ind)
  x<-x[!na.ind]
  n<-length(x)
  if((n / (n+na.n)) < 0.9)
    return(NA)
  y<-sum(x>=h_thresh)
  return(y)
}

noncover_n<-function(x, h_thresh, na.rm){
  na.ind<-is.na(x)
  na.n<-sum(na.ind)
  x<-x[!na.ind]
  n<-length(x)
  if((n / (n+na.n)) < 0.9)
    return(NA)
  y<-sum(x<h_thresh)
  return(y)
}

cover_resid_calc<-function(cover, tch){
  cover_resid<-cover - (1 + exp(12.431) * tch^-4.061)^-1
  return(cover_resid)
}


#' Calculate basal area from top canopy height and the residual of the relationship
#' between tch and canopy cover proportion
#'
#' @description Functions for calculating AGB taken from Jucker et al., 2017
#' "A regional model for estimating the aboveground carbon density
#' of Borneo's tropical forests from airborne laser scanning"

#' @param x A raster image
#' @param outshape The name of the shapefile if saving
#' @param gdal_format The format of the shapefile to use; by default ESRI shapefile format
#' @param pypath The path to gdal_polygonize.py on the system
#' @param readpoly Logical, encoding whether or not the output should be returned or not. If
#' not NULL is returned.
#' @param quiet Logical indicating whether or not a verbose output should be creaeted.
#' @return A spatialPolygonsDataframe
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-09-13
ba <- function(tch, cover_resid){
  1.287 * tch^0.987 * (1 + 1.983 * cover_resid)
}


#' Calculate wood density from top canopy height
#'
#' @description Functions for calculating AGB taken from Jucker et al., 2017
#' "A regional model for estimating the aboveground carbon density
#' of Borneo's tropical forests from airborne laser scanning"

#' @param tch Top canopy height as a numeric vector
#' @return A numeric vector of wood densities
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-09-13
wd <-function(tch){
  0.385 * tch^0.097
}


#' Estimate above ground biomass
#'
#' @description Functions for calculating AGB taken from Jucker et al., 2017
#' "A regional model for estimating the aboveground carbon density
#' of Borneo's tropical forests from airborne laser scanning"

#' @param tch Top canopy height as a numeric vector
#' @param cover_residual Canopy cover proportion residual as a numeric vector
#' @return above ground biomass in tonnes per ha as a numberic vector
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-09-13
est_AGB<-function(tch, cover_resid) {
  agb <- 0.567 * tch^0.554 * ba(tch, cover_resid)^1.081 * wd(tch)^0.186
  return(agb)
  #1.44*tch^0.63 * (6.44*tch)^0.85 * (0.489 + 0.0022*tch)^0.24
}


#' Find the nearest negative integer
#'
#' @param x A numeric vector
#' @return A vector of negative integers
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-09-13
nearest_neg_round<-function(x)
  2*round((x+1)/2)-1


#' Create a stack of tch and canopy cover:
#'
#' @description Create a stack of tch and canopy cover:
#'
#' @param chm A raster image of canopy height values
#' @param cover_h_thresh A numeric value indicating the height at which canopy cover proportion
#' should be calculated
#' @param scale The scale over which agb should be calculated in meters
#' @param verbose logical; indicates whether or not progress should be printed to the console
#' @return A raster with pixel values denoting above ground biomass in tonnes per ha. Note
#' that these values should not be added in aggregation but averaged instead.
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 18-03-27
#'
tch_stacker <- function(chm, cover_h_thresh, scale, verbose = TRUE) {

  n_pix<-scale / res(chm)[1]
  n_pix<-nearest_neg_round(n_pix)
  w<-matrix(1, n_pix, n_pix)
  progress<-""
  if(verbose) progress<-"text"
  #lid_cover<-focal(lid_chm, w = matrix(1, 301, 301), fun = function(x) cover(x, h_thresh = 20))
  if(verbose) cat("Calculating cover: \n")
  cover_n <-
    aggregate(
      chm,
      fact = n_pix,
      fun = function(x, na.rm)
        cover_n(x, h_thresh = cover_h_thresh),
      progress = progress
    )
  if(verbose) cat("\nCalculating openness: \n")
  noncover_n <-
    aggregate(
      chm,
      fact = n_pix,
      fun = function(x, na.rm)
        noncover_n(x, h_thresh = cover_h_thresh),
      progress = progress
    )
  if(verbose) cat("\nCalculating top canopy height: \n")
  tch <- aggregate(chm, fact = n_pix, fun = mean, progress = progress)

  tch_stack<-stack(tch, cover_n, noncover_n)
  names(tch_stack)<-c("tch", "cover", "gap")
  return(tch_stack)
}




#' Calculate a raster of above ground biomass
#'
#' @description Functions for calculating AGB taken from Jucker et al., 2017
#' "A regional model for estimating the aboveground carbon density
#' of Borneo's tropical forests from airborne laser scanning"
#'
#' @param chm A raster image of canopy height values
#' @param cover_h_thresh A numeric value indicating the height at which canopy cover proportion
#' should be calculated
#' @param scale The scale over which agb should be calculated in meters
#' @param verbose logical; indicates whether or not progress should be printed to the console
#' @return A raster with pixel values denoting above ground biomass in tonnes per ha. Note
#' that these values should not be added in aggregation but averaged instead.
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 17-09-13
raster_agb <- function(chm, cover_h_thresh, scale, est_tch_gf = FALSE, verbose = TRUE) {

  n_pix<-scale / res(chm)[1]
  n_pix<-nearest_neg_round(n_pix)
  w<-matrix(1, n_pix, n_pix)
  progress<-""
  if(verbose) progress<-"text"
  #lid_cover<-focal(lid_chm, w = matrix(1, 301, 301), fun = function(x) cover(x, h_thresh = 20))
  if(verbose) cat("Calculating cover: \n")
  cover_n <-
    aggregate(
      chm,
      fact = n_pix,
      fun = function(x, na.rm)
        cover_n(x, h_thresh = cover_h_thresh),
      progress = progress
    )
  if(verbose) cat("\nCalculating openness: \n")
  noncover_n <-
    aggregate(
      chm,
      fact = n_pix,
      fun = function(x, na.rm)
        noncover_n(x, h_thresh = cover_h_thresh),
      progress = progress
    )
  if(verbose) cat("\nCalculating top canopy height: \n")
  tch <- aggregate(chm, fact = n_pix, fun = mean, progress = progress)

  if(verbose) cat("\nEstimating AGB: \n")
  tf = cbind(cover_n = values(cover_n),
             cover_notn = values(noncover_n))

  biomass_df <-
    data.frame(
      tch = values(tch),
      n = values(cover_n),
      notn = values(noncover_n)
    ) %>%
    mutate(cover = n / (n + notn))

  if(est_tch_gf){
    ggplot(biomass_df, aes(x = tch, y = cover)) +
      geom_point() +
      geom_smooth(method = "glm",
                  method.args = list(family = "binomial"))

    fm.cover1 <<- glm(tf ~ tch, data = biomass_df, family = binomial)
    #summary(fm.cover1)
    biomass_df$cover_resid <- NA
    biomass_df$cover_resid[!is.na(biomass_df$n)] <-
      residuals(fm.cover1, type = "response")
  }
  else{
    biomass_df$cover_resid <-cover_resid_calc(biomass_df$cover, biomass$tch)
  }

  agb <- cover_n
  agb[] <- est_AGB(biomass_df$tch, biomass_df$cover_resid)

  return(agb)
}
