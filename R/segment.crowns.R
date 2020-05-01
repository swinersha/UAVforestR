#' Segments crowns for a set of maxima coordinates
#'
#' @description Takes the coordSeeds object created by detect.maxima as well as the image and
#' sobel images and segments out the crowns.
#'
#' @param x The image encoded as a matrix (this is usally the blurred version of the
#' original image)
#' @param x.Sobel The Sobel edge image encoded as a matrix.
#' @param coordSeeds A coordSeeds object creaeted by detect.maxima()
#' @param THRESHSeed The crown minimum height as a proportion (which is relative to
#' maximum crown height)
#' @param THRESHCrown The crown minimum height as a proportion (which is relative to
#' mean crown height, which is adjusted on the fly during processing)
#' @param SOBELstr The sobel stength used to detect edges; this is a little cryptic
#' to interpret. It is set relative to the height of the tree and essentially
#' encodes the edge strength needed to define an edge. Larger values of SOBELstr
#' will make the edge detect less sensitive. Values between 0 and 5 are usually a good
#' starting point.
#' @param scale The image resolution in meters
#' @param tau The percentile to use in calculating the tree maximum size.
#' @return A list containing: 1. The pixel locaitons within the matrix by tree id;
#' 2. The pixel locations of the tree maximum, along with tree heights, the proportion of the
#' crown edge detected by each boundary condition, and the average sobel strength of the
#' tree boundary.
#' @export
#' @author Tom Swinfield
#' @details Created 17-08-17

segment.crowns <- function(x, x.Sobel, coordSeeds, THRESHSeed, THRESHCrown, SOBELstr, scale, lut, tau) {

  Check <- x #
  Check[, ] <- 0 # Sets check to 0.

  length.total<-ncell(x)

  coordSeeds<-coordSeeds[order(coordSeeds[,'height'], decreasing = TRUE),]
  treeHeights_mean<-rep(0, length = nrow(coordSeeds))
  treeHeights_max<-rep(0, length = nrow(coordSeeds))
  sobel_mean<-vector('numeric', length=nrow(coordSeeds))
  boundcount<-matrix(0, nrow=nrow(coordSeeds), ncol=6, dimnames=list(NULL, c('allombound', 'crownbound', 'hbound_max', 'hbound_mean', 'sobelbound', 'tallerbound')))

  coordSeeds<-coordSeeds[,c('row', 'col')]

  crown.ind<-1 # A marker to keep track of the active pixel
  crown.total<-crown.ind  # A marker to keep track of the number of pixels still left to work on.

  coordCrown<-matrix(NA, nrow=length.total, ncol=2, dimnames=list(NULL, c('row', 'col')))
  # sobel_mean<-vector('numeric', length=length.total)
  # boundcount<-matrix(0, nrow=length.total, ncol=5, dimnames=list(NULL, c('hbound', 'sobelbound', 'crownbound', 'allombound', 'tallerbound')))
  tree.ind<-vector(length.total, mode='numeric')

  Ntrees<-nrow(coordSeeds)

  for (i in 1:Ntrees)
  {
    newCrown <-
      matrix(coordSeeds[i, ],
             ncol = 2,
             dimnames = list(NULL, c('row', 'col'))) # The pixel for the new crown seed is marked as newCrown, to initiate the next iteration.

    if(!Check[newCrown]){
      coordCrown[crown.ind, ] <- newCrown

    # Check[r, k + 1]

    crown.indSeed <-
      crown.ind # The index of the first pixel of the crown.
    coordSeed <- newCrown # The treeseed for the current crown.

    rvSeed <- x[coordSeed] # Extracts the tree height.
    rvSobel <-
      # rvSeed * (1 - THRESHCrown) * SOBELstr # This multiplier changes the sensitivity of the sobel segment;
      SOBELstr # This multiplier changes the sensitivity of the sobel segment;
    # it might be better if this is set to THRESHCrown

    Check[newCrown] <- i # and is checked off.
    rvCrown <-
      mean(x[newCrown], na.rm = T) # Extracts the mean tree height.
    #crownrad.THRESH<-(htod(rvCrown, lut, tau=90)/2)/res(imagery)[1] # The crown radius in pixels
    crownrad.THRESH <-
      (htod(x = rvCrown, lut = lut, tau = tau) / 2) / scale # The crown radius in pixels
    boundMat <- boundMat_set(coordSeed[1], coordSeed[2], rvCrown, scale = scale, lut = lut, tau = tau)

    crown.total <-
      crown.ind  # A marker to keep track of the number of pixels still left to work on.

    it.chk <- 0
    while (crown.ind <= crown.total)
    {
      it.chk <- it.chk + 1
      r = as.numeric(coordCrown[crown.ind, 1])
      k = as.numeric(coordCrown[crown.ind, 2]) # Extracts the tree location.
      #if (r != 1 & r != dim(x)[1] & k != 1 & k != dim(x)[2])
      # { # So longs as the pixel is not right on the boundary; !!! this might need to be changed depending on window size.

      #rvCrown <- mean(x[coordCrown[,1]+nrow(x)*(coordCrown[,2]-1)], na.rm = T) # Extracts the mean tree height.
      crownHeights <-
        x[matrix(coordCrown[crown.indSeed:crown.ind, ],
                    ncol = 2,
                    dimnames = list(NULL, c('row', 'col')))]
      rvCrown <- mean(crownHeights) # Extracts the mean tree height.
      rvCrown_sd <- sd(crownHeights)
      # You could make this more efficient using weighted means but whatever; it's way less readable.

      # Makes a matrix of the coordinates and chm values for the adjacent cells (surrounding the focal pixel).
      filData <- matrix(4, 5, data = 0)
      filData[1, 1] <- r - 1
      filData[1, 2] <- k
      filData[1, 3] <- x[r - 1, k]
      filData[1, 4] <- Check[r - 1, k]
      filData[1, 5] <- x.Sobel[r - 1, k]
      filData[2, 1] <- r
      filData[2, 2] <- k - 1
      filData[2, 3] <- x[r, k - 1]
      filData[2, 4] <- Check[r, k - 1]
      filData[2, 5] <- x.Sobel[r, k - 1]
      filData[3, 1] <- r
      filData[3, 2] <- k + 1
      filData[3, 3] <- x[r, k + 1]
      filData[3, 4] <- Check[r, k + 1]
      filData[3, 5] <- x.Sobel[r, k + 1]
      filData[4, 1] <- r + 1
      filData[4, 2] <- k
      filData[4, 3] <- x[r + 1, k]
      filData[4, 4] <- Check[r + 1, k]
      filData[4, 5] <- x.Sobel[r + 1, k]


      # Calculates distance of pixels from the focal pixel
      fil.dists <-
        as.matrix(dist(rbind(coordSeed, filData[, 1:2])))[1, -1]

      # Edits the boundary conditions matrix
      boundMat <-
        boundMat_edit(
          boundMat,
          filData,
          fil.dists,
          i,
          rvSeed,
          rvSobel,
          THRESHCrown,
          THRESHSeed,
          crownrad.THRESH
        )
      #boundMat_view(boundMat_edit(boundMat, filData, fil.dists))

      # Calculates the crown confidence envelope
      Crown_upper <- rvCrown + rvCrown_sd
      Crown_lower <- rvCrown - rvCrown_sd

      # Checks which (if any) of the values are greater than the max or mean tree heights adjusted by the thresholds.
      # and less than 5% taller than the max tree height.
      # and the euclidean distance from the maximum crown radius is less than the specified DIST.
      GFIL <- (
        filData[, 3] > (rvSeed * THRESHSeed) &
          (filData[, 3] > (rvCrown * THRESHCrown)) &
          (filData[, 3] <= (rvSeed + (rvSeed * 0.05))) &
          (fil.dists < crownrad.THRESH) &
          (filData[, 4] == 0) &
          (filData[, 5] < rvSobel)
      )
      newCrown <-
        matrix(filData[GFIL, 1:2],
               ncol = 2,
               dimnames = list(NULL, c('row', 'col'))) # Subsets filData by the decision tree output.

      nNew <- sum(GFIL) # The number of new crown pixels
      # Adds the pixels that pass the test to the crown and checks them off as done.
      if (nNew > 0)
        # Does this need to be 3?
      {
        Check[newCrown] <- i
        coordCrown[(crown.total + 1):(crown.total + nNew), ] <- newCrown
        crown.total <-
          crown.total + nNew # The total size of the crown data we're working on.
      }
      crown.ind <- crown.ind + 1 # increment focal pixel
    }
    } # End of crown while loop.
    tree.ind[crown.indSeed:(crown.ind - 1)] <- i
    #treeHeights[ind] <- max(crownHeights) # Extracts the max tree height.
    treeHeights_mean[i] <-
      mean(crownHeights) # Extracts the mean tree height.
    treeHeights_max[i] <-
      max(crownHeights) # Extracts the mean tree height.
    boundcount[i,]<-boundMat_summarise(boundMat)
    sobel_mean[i]<-boundMat_sobel_mean(boundMat, x.Sobel)
  }

  # Extracts only the filled parts of the data objects:
  coordCrown<-coordCrown[1:crown.total,]
  coordCrown <- cbind(coordCrown, id = tree.ind[1:crown.total]) # Extracts the tree crowns, which are already in order.

  boundcount<-boundcount[1:Ntrees,]
  boundcount<-prop.table(boundcount, margin = 1) # convert to proportions
  sobel_mean<-sobel_mean[1:Ntrees]
  coordSeeds <- cbind(coordSeeds,
                      treeHeights_max,
                      treeHeights_mean,
                      boundcount,
                      sobel_mean = sobel_mean
                      )
  coordSeeds<-as.data.frame(coordSeeds)

  # Create output object
  coordCrown<-list(coordCrown, coordSeeds)

  #figure out what needs to be returned
  return(coordCrown)
}
