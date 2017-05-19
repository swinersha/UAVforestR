##########################################################
# A modification of ITCsegment package by Michele Dalponte

# Changes made by Tom Swinfield for the RSPB under funding from the Cambridge Conservation Initiative
# Collaborative Fund
##########################################################

# This is a reverse engineering of ITC segment topdown to first find all the local maxima and
# then to work through these from tallest to shortest using the fast algorithm.

# Changes made:

# 160607: - Enabled the search window to size to be scalable.
#         - Tidied up some of the code.
#         - Made the maxima search tool work from tallest to smallest.

# 160711: - Changed the distance of pixels to the local maximum a raial distance rather than straight line distance.
#         - Scale the maximum distance according to the allometry.
#           The programme should work through the trees from tallest to shortest and identify the crown area on the fly; 
#           This would over come problems associated with climbing up the canopy on to adjacent crowns.

# 161010: - WARNING: if gobble is on; the mean heights will be incorrect.
#           The gobble function will have to recalculate the means, probably by weighting them by crown size if they are
#           to be accurate.


# You need to make sure you have GDAL installed and that the file path to gdal_polygonize.py is set correctly

# You need to make it easy to put in an allometric equation and then scale everything relative to that.


#if (r != 1 & r != dim(Gnew)[1] & k != 1 & k != dim(Gnew)[2]) 
# So longs as the pixel is not right on the boundary; !!! 
# this will need to be changed as the window size is modified.

# Is the gaussian blur necessary?

#THRESHSeed=0.45 # Proportional height relative to tree maximum.
#THRESHCrown=0.55 # Proportional height relative to tree mean.
#DIST=8 # Distance in pixels from tree maximum.
#specT=2 # Minimum height in m for tree or crown.

# x is the height
# htod_lut is a lookup table containing the allometric dog leg coeficients
# tau should be entered as the percentile.
htod_dc<-function(x, htor_lut, tau) 
{
  rad<-x
  a<-htor_lut[tau,'a']
  b<-htor_lut[tau,'b']
  aa<-htor_lut[tau,'aa']
  bb<-htor_lut[tau,'bb']
  rad[x<exp(3)]<-(exp(aa) * x[x<exp(3)]^bb )
  rad[x>=exp(3)]<-(exp(a) * x[x>=exp(3)]^b)
  return(rad*2)
}

htod_hrf<-function(x, htor_lut, tau) 
{
  rad<-x
  a<-htor_lut[tau,'aa']
  b<-htor_lut[tau,'bb']
  rad<-(exp(a) * x^b)
  return(rad*2)
}

htod_lookup<-function(x, lut, tau) 
{
  rad<-x
  a<-lut[tau,'a']
  b<-lut[tau,'b']
  rad<-(exp(a) * x^b) 
  return(rad*2)
}

# Calculates the window size when given the relationship between tree height and window size as a fitted model (fm.win)
# and z the tree heights:
allom.dist<-function(z, scale)
{
  diam<-htod(z)
  diam<-(diam/scale)
  win.size<-2*round((diam+1)/2)-1
  win.size[win.size<3]<-3
  names(win.size)<-NULL
  return(win.size)
}
# A function to calculate the number of pixels distance to the edge
edge.dist<-function(z, scale)
{
  radius<-htod(z)/2
  radius<-radius/scale
  pixels<-ceiling(radius)
  names(pixels)<-NULL
  return(pixels)
}

#imagery=uav_chm_blur
#im_sobel=uav_chm_sobel
#THRESHSeed=0.7; THRESHCrown=0.7; specT=2; lm.searchwin=NULL; SOBELstr=2; can_text=NULL; cansize=20
#htod=function(x) htod_lookup(x, lut=lut, tau=90)

# Extracts just the locations of the tree maxima:
itcIMG_fast<-function (imagery, im_sobel=NULL, THRESHSeed, 
  THRESHCrown, htod, specT, SOBELstr, lm.searchwin=NULL, gobble='off', cantext='off', cansize=NULL)
{
  # !!! You need to make some appropriate checks of the arguments here !!!
  if(cantext=='on' & is.null(cansize))
    stop('In itcIMG_fast cansize must not be NULL if cantext is on; please enter numberic value')
  
  cat('Extracting matrix\n')
  # Extracts the image data as a matrix:
  Max <- matrix(dim(imagery)[2], dim(imagery)[1], data = imagery[,], byrow = FALSE)
  Sobel <- matrix(dim(imagery)[2], dim(imagery)[1], data = im_sobel[,], byrow = FALSE)
  cat('Preparing data structures\n')
  # Flips the image the right way again:
  Max <- Max[1:dim(imagery)[2], dim(imagery)[1]:1]
  Sobel <- Sobel[1:dim(imagery)[2], dim(imagery)[1]:1]
  Gnew <- Max     # Copies the max matrix.
  Max[, ] <- 0    # Sets max to 0s. This will be used later for storing ...
  Index <- Max    # Copies max again.
  Index[, ] <- 0  # Sets index to 0s.
  Gnew[is.na(Gnew)] <- 0 # Sets any nas to 0s.
  Gnew[Gnew < specT] <- 0 # Any values beneath the minimum height are set to 0s.
  index = 1       # Initiates the index
  II <- which(Gnew != 0, arr.ind = T) # Extracts the locations of pixels which are bigger than the min tree height.
  dim(II)
  
  #fm.win<-win.lmfun(win.min, win.max, z.min, z.max) # generates the relationship between height and search window size.
  # Extracts only the pixels that are sufficiently far from the image edge for the search window. In each direction.
  # The search window is selected according to the height of the tree:
  z<-Gnew[II[,1]+nrow(Gnew)*(II[,2]-1)] # extracts the tree heights from the matrix
  #z.pix<-z/res(imagery)[1] # converts heights in m to pixels
  #WinSize<-allom.dist(z.pix) # Finds the window size for each tree pixel in m
  WinSize<-allom.dist(z, scale=res(imagery)[1]) #, lut, tau=10)
  pix2edge<-edge.dist(z, scale=res(imagery)[1])#, lut, tau=90) # Halfs the window size for subsetting by pixels far enough from the image edge
  
  # Extracts only the pixels far enough from the image edge.
  II.ind<-II[, 1] > pix2edge & 
    II[, 1] < (nrow(Gnew) - pix2edge) &
    II[,2] > pix2edge &
    II[,2] < (ncol(Gnew) - pix2edge)
  II<-II[II.ind,]
  WinSize<-WinSize[II.ind]
  z<-z[II.ind]
  dim(II)
  
  cat('Ordering pixels\n')
  # reorder from greatest to least according to z value (i.e. big trees first):
  z.order<-order(z, decreasing = TRUE)
  II<-II[z.order,]
  WinSize<-WinSize[z.order]
  
  # Extracts the chm values:
  Cb <- imagery
  Mb <- imagery
  Cb[] <- as.numeric(Gnew[1:dim(Gnew)[1], dim(Gnew)[2]:1], 
                     byrow = TRUE)
  Mb[] <- as.numeric(Max[1:dim(Max)[1], dim(Max)[2]:1], 
                     byrow = TRUE)
  
  Crowns <- Index # Assigns the crown sequence within the raster to Crowns
  OldCrowns <- Crowns # Copies crowns
  Check <- OldCrowns # Copies again.
  Check[, ] <- 0 # Sets check to 0.
  
  ind<-1 # This is the tree index and should just increment from 1 (largest) to n (the smallest tree).
  crown.ind<-1 # A marker to keep track of the active pixel
  crown.total<-crown.ind  # A marker to keep track of the number of pixels still left to work on.
  length.total<-length(Cb)
  length.II<-nrow(II)
  coordCrown<-matrix(NA, nrow=length.total, ncol=2, dimnames=list(NULL, c('row', 'col')))
  coordSeeds<-matrix(NA, nrow=length.total, ncol=2, dimnames=list(NULL, c('row', 'col')))
  sobel_mean<-vector('numeric', length=length.total)
  treeHeights<-vector('numeric', length=length.total)
  boundcount<-matrix(0, nrow=length.total, ncol=5, dimnames=list(NULL, c('hbound', 'sobelbound', 'crownbound', 'allombound', 'tallerbound')))
  tree.ind<-vector(length.total, mode='numeric')
  
  # Works through each pixel one by one:
  cat('Finding trees\n')
  #indexII<-1
  #indexII<-(1:nrow(II)[1])[1]
  for (indexII in 1:dim(II)[1]) 
  {
    r = as.numeric(II[indexII, 1]) # Extracts the row pos.
    k = as.numeric(II[indexII, 2]) # Extracts the column pos.
    # Sets the search window size for the local maxima:
    if(class(lm.searchwin)=='NULL') 
      searchWinSize<-WinSize[indexII] # a search window scaled by max tree height.
    else
      searchWinSize<-lm.searchwin # a fixed size search window.
    
    hc.sWS<-ceiling(searchWinSize/2) # half the search window size rounded to floor
    hf.sWS<-floor(searchWinSize/2) # half the search window size rounded to ceiling
    FIL <- matrix(searchWinSize, searchWinSize, data = NA) 
    # Extracts the search window for the pixel:
    FIL <- Gnew[(r - hf.sWS):(r + hf.sWS), 
                (k - hf.sWS):(k + hf.sWS)]
    # Extracts the window from Max indicating whether a tree has already being designated or not:
    Max.chk<-Max[(r - hf.sWS):(r + hf.sWS),
                 (k - hf.sWS):(k + hf.sWS)]
    
    # If the focal pixel has the greatest value in the window & there is no tree already assigned in the output matrix within the window & the max value is no 0...
    # because the order is from tallest to shortest, large trees will always suppress the designation of small trees.
    if (FIL[hc.sWS, hc.sWS] == max(FIL, na.rm = T) &
        max(Max.chk, na.rm = T) == 0 &
        max(FIL, na.rm = T) != 0) 
    {
      Max[r, k] <- 1 # A logical assignment of tallest tree
      coordSeeds[index, ] <- c(r,k) # Assigns the sequence in which the trees were found.
      treeHeights[index]<-Gnew[r,k] # Assigns the height of the trees
      index <- index + 1 # increments
    }
  }
  Ntrees <- index -1 # Number of trees encountered.
  cat(Ntrees, 'trees detected; creating spatial objects\n')  
  coordSeeds<-coordSeeds[1:Ntrees,] # Extracts the tree crowns, which are already in order.
  treeHeights<-treeHeights[1:Ntrees] # Extracts the tree crowns, which are already in order.
  treeHeights_mean<-treeHeights # Copies treeHeights so that the mean crown values can be saved.

  cat('Segmenting crowns\n')  
  #ind<-1
  for(ind in 1:Ntrees)
  {  
      newCrown<-matrix(coordSeeds[ind,], ncol=2, dimnames=list(NULL, c('row', 'col'))) # The pixel for the new crown seed is marked as newCrown, to initiate the next iteration.  
      coordCrown[crown.ind,]<-newCrown
      crown.indSeed<-crown.ind # The index of the first pixel of the crown.
      coordSeed<-newCrown # The treeseed for the current crown.
      Check[newCrown]<-ind # and is checked off.
      rvCrown<-mean(Gnew[newCrown], na.rm = T) # Extracts the mean tree height.
      #crownrad.THRESH<-(htod(rvCrown, lut, tau=90)/2)/res(imagery)[1] # The crown radius in pixels
      crownrad.THRESH<-(htod(x=rvCrown)/2)/res(imagery)[1] # The crown radius in pixels
      
      crown.total<-crown.ind  # A marker to keep track of the number of pixels still left to work on.
      
      it.chk<-0
      while (crown.ind<=crown.total)
      {
        it.chk<-it.chk+1
        r = as.numeric(coordCrown[crown.ind, 1])
        k = as.numeric(coordCrown[crown.ind, 2]) # Extracts the tree location.
        #if (r != 1 & r != dim(Gnew)[1] & k != 1 & k != dim(Gnew)[2]) 
        # { # So longs as the pixel is not right on the boundary; !!! this might need to be changed depending on window size.
        
        rvSeed <- Gnew[coordSeed] # Extracts the tree height.
        rvSobel <- rvSeed * (1-THRESHSeed) * SOBELstr # This multiplier changes the sensitivity of the sobel segment; 
        # it might be better if this is set to THRESHCrown
        #rvCrown <- mean(Gnew[coordCrown[,1]+nrow(Gnew)*(coordCrown[,2]-1)], na.rm = T) # Extracts the mean tree height.
        crownHeights<-Gnew[matrix(coordCrown[crown.indSeed:crown.ind,], ncol=2, dimnames=list(NULL, c('row', 'col')))]
        rvCrown <- mean(crownHeights) # Extracts the mean tree height.
        rvCrown_sd<- sd(crownHeights)
        # You could make this more efficient using weighted means but whatever; it's way less readable.
        
        # Makes a matrix of the coordinates and chm values for the adjacent cells (surrounding the focal pixel).
        filData <- matrix(4, 5, data = 0)
        filData[1, 1] <- r - 1
        filData[1, 2] <- k
        filData[1, 3] <- Gnew[r - 1, k]
        filData[1, 4] <- Check[r - 1, k]
        filData[1, 5] <- Sobel[r - 1, k]
        filData[2, 1] <- r
        filData[2, 2] <- k - 1
        filData[2, 3] <- Gnew[r, k - 1]
        filData[2, 4] <- Check[r, k - 1]
        filData[2, 5] <- Sobel[r, k - 1]
        filData[3, 1] <- r
        filData[3, 2] <- k + 1
        filData[3, 3] <- Gnew[r, k + 1]
        filData[3, 4] <- Check[r, k + 1]
        filData[3, 5] <- Sobel[r, k + 1]
        filData[4, 1] <- r + 1
        filData[4, 2] <- k
        filData[4, 3] <- Gnew[r + 1, k]
        filData[4, 4] <- Check[r + 1, k]
        filData[4, 5] <- Sobel[r + 1, k]
        
        filData.CHK<<-filData
        # Calculates distance of pixels from the focal pixel
        fil.dists<-as.matrix(dist(rbind(coordSeed, filData[,1:2])))[1,-1]
        # Calculates the crown confidence envelope
        Crown_upper <- rvCrown + rvCrown_sd
        Crown_lower <- rvCrown - rvCrown_sd
        
        # Checks which (if any) of the values are greater than the max or mean tree heights adjusted by the thresholds.
        # and less than 5% taller than the max tree height.
        # and the euclidean distance from the maximum crown radius is less than the specified DIST.
        GFIL <- (filData[, 3] > (rvSeed * THRESHSeed) & 
                   (filData[, 3] > (rvCrown * THRESHCrown)) & 
                   (filData[, 3] <= (rvSeed + (rvSeed * 0.05))) &
                   (fil.dists<crownrad.THRESH) &
                   (filData[,4] == 0) &
                   (filData[,5] < rvSobel))
        newCrown <- matrix(filData[GFIL, 1:2],ncol=2, dimnames=list(NULL, c('row', 'col'))) # Subsets filData by the decision tree output.
        
        # Storing some information about the reasons for segmentation according to the parameters:
        HB<-filData[, 3] <= (rvSeed * THRESHSeed) # pixel less than max height threshold
        MHB<-filData[, 3] <= (rvCrown * THRESHCrown) # pixel less than mean height threshold
        SEB<-filData[, 5] >= rvSobel # Pixel greater than edge threshold
        CRB<-!(filData[,4] %in% c(0,ind)) # pixel in another tree crown
        ADB<-fil.dists >= crownrad.THRESH # pixel greater than allometric threshold
        IHB<-filData[, 3] > (rvSeed + (rvSeed * 0.05)) # pixel taller than tree maximum
        
        # Workds out the number of boundary conditions for each:
        nhb<-sum(HB | MHB) # Height boundary
        nsb<-sum(SEB & !(HB | MHB)) # Sobel boundary
        ncb<-sum(CRB & !(HB | MHB| SEB)) # Crown boundary
        nab<-sum(ADB & !(HB | MHB | SEB | CRB)) # Allometric boundary
        ntb<-sum(IHB & !(HB | MHB | SEB | CRB | ADB)) # Taller tree boundary
        
        if(any(SEB | HB | MHB))
        {
          sobel_weights<-c(sum(boundcount[ind,3:5]), rep(1, sum(nhb+nsb)))
          sobel_mean[ind]<-weighted.mean(c(sobel_mean[ind], filData[(SEB | HB | MHB), 5]), w=sobel_weights)
        }
        
        boundcount[ind,1]<- boundcount[ind,1] + nhb
        boundcount[ind,2]<- boundcount[ind,2] + nsb
        boundcount[ind,3]<- boundcount[ind,3] + ncb
        boundcount[ind,4]<- boundcount[ind,4] + nab
        boundcount[ind,5]<- boundcount[ind,5] + ntb
        
        nNew<-sum(GFIL) # The number of new crown pixels
        # Adds the pixels that pass the test to the crown and checks them off as done.
        if (nNew > 0) # Does this need to be 3?
        {
          Check[newCrown]<-ind
          coordCrown[(crown.total+1):(crown.total+nNew),]<-newCrown
          crown.total<-crown.total+nNew # The total size of the crown data we're working on.
        }
        crown.ind<-crown.ind+1 # increment focal pixel
      } # End of crown while loop.
      tree.ind[crown.indSeed:(crown.ind-1)]<-ind
      #treeHeights[ind] <- max(crownHeights) # Extracts the max tree height.
      treeHeights_mean[ind] <- mean(crownHeights) # Extracts the mean tree height.
#    } # If end
#    else {continue<-0}
  } # End of tree for loop.
  
  # Extracts only the filled parts of the data objects:
  tree.ind<-tree.ind[1:crown.total]
  coordCrown<-coordCrown[1:crown.total,] 
  boundcount<-boundcount[1:Ntrees,]
  sobel_mean<-sobel_mean[1:Ntrees]
  
  # Convert boundary condition counts to proportions:
  boundcount<-prop.table(boundcount[,1:5], 1)
  # !!!!! This needs to be fixed !!!!!
  
  #plot(y=sobel_mean,x=treeHeights)
  #plot(y=sobel_mean/(4*treeHeights),x=rowSums(boundcount[,3:4]))
  #hist(sobel_mean/(4*treeHeights), breaks=30)
  #plot(coordCrown, col=tree.ind, pch=16)
  
  if(gobble=='on')
  {
    cat('Gobbling\n')
  # This requires a few objects:
    # 1. tree.ind, which is a vector of the identities of the crowns by pixel, 
    # which should match up with a which(Crowns, arr.ind=TRUE) object for the
    # whole Crowns matrix. - I guess this should be easy for you to make from the Crowns
    # matrix.
    # 2. coordSeeds, which is a which(Index, arr.ind=TRUE) object for each of the
    # maxima from tallest to shortest.
    # 3. treeHeights, which is a vector of the tree heights for each of the maxima
    # 4. the function htod, which provides the look up table for allometry coefficients.
    # 5. the parameters: 
      # hthresh_gob=20 - the min height threshold for gobbling
      # hpropmin_gob=1 - the min height for a bigger tree to combine with
      # hpropmax_gob=1.5 - the max height for a bigger tree to combine with
      # alpharad_gob=2 - a expansion coefficient for the 90% radius for the tree height * hpropmax_gob
    
    
    treeAreas<-as.numeric(table(tree.ind)) * res(imagery)[1]^2
    treeRadii<-sqrt(treeAreas/pi)
    treeRadii_lower<-htod(treeHeights, lut, tau=20)/2
    treeRadii_upper<-htod(treeHeights, lut, tau=90)/2
    #plot(y=treeRadii, x=treeRadii_lower); abline(0,1)
  
    #plot(treeRadii~treeHeights)
    #points(htod(treeHeights, lut, tau=20)/2~treeHeights, type='p', col='blue')
    #points(htod(treeHeights, lut, tau=50)/2~treeHeights, type='p', col='red')
    #points(htod(treeHeights, lut, tau=90)/2~treeHeights, type='p', col='purple')
  
  #  Finds the distance between all tree maxima
    treedist<-as.matrix(dist(coordSeeds, upper=TRUE))
    diag(treedist)<-NA
  
  # checks which trees are smaller and bigger than the 10% percentile respectively:
  
    treeID_toosmall<-which(treeRadii_lower>treeRadii & treeHeights>hthresh_gob)
    #treeID_toosmall<-which(treeRadii==sqrt((1 * res(imagery)[1]^2)/pi)) # 1 pixel
    treeID_rightsize<-which(treeRadii_lower<=treeRadii & treeHeights>hthresh_gob)
  
    points(y=treeRadii[treeID_toosmall], x=treeHeights[treeID_toosmall], col='green')
  
    treedist.sub<-treedist[treeID_toosmall,treeID_rightsize]
    treeHeights_toosmall<-treeHeights[treeID_toosmall]
    treeHeights_rightsize<-treeHeights[treeID_rightsize]
    treeAreas_toosmall<-treeAreas[treeID_toosmall]
    treeAreas_rightsize<-treeAreas[treeID_rightsize]
  
    for(i in length(treeID_toosmall):1)
    {
    # finds the trees within a suitable radius distance:
      min_height<-treeHeights_toosmall[i]*hpropmin_gob
      max_height<-treeHeights_toosmall[i]*hpropmax_gob
      max_dist<-(htod90(max_height)/2)*alpharad_gob
      match_treeID<-treedist.sub[i,]<=max_dist & 
        treeHeights_rightsize>min_height &
        treeHeights_rightsize<=max_height
      match_tree<-treeID_rightsize[match_treeID]
      if(length(match_tree)>=1)
      {
      # Find the tree with the biggest attraction:
        G<-(treeAreas_rightsize[match_treeID]*treeAreas_toosmall[i])/(treedist.sub[i,match_treeID]^2) # calculates the gravitational force.
      
        gobbleareas<-treeAreas_rightsize[match_treeID]
        gobbleHeights<-treeHeights_rightsize[match_treeID]
        gobbleareas_upper<-pi*(htod(gobbleHeights, lut, tau=90)/2)^2
        combarea_test<-(gobbleareas + treeAreas_toosmall[i])<=gobbleareas_upper # checks if the combined areas woule exceed the 99% percentile
        match_tree<-match_tree[which(combarea_test)]   # any that are, are excluded.
        G<-G[which(combarea_test)]                     # and also excluded from here.
        gobbleareas<-gobbleareas[which(combarea_test)] # and also excluded from here.
        gobbleareas<-gobbleareas[which.max(G)] # extracts the area of the gobbling tree.
        gobbleID<-match_tree[which.max(G)] # finds the single tree with the biggest gravitational force.
      
        if(length(gobbleID)==1)
        {
          if(length(gobbleID)>1)
            cat(i, 'Woah too many gobblers\n')
          treeID<-treeID_toosmall[i]
          tree.ind[tree.ind==treeID]<-gobbleID
        
          gobbleareas_upper<-gobbleareas_upper[which(combarea_test)]
          gobbleareas_upper<-gobbleareas_upper[which.max(G)]
          newarea<-sum(tree.ind==gobbleID) * res(imagery)[1]^2

          if(newarea>gobbleareas_upper)
            cat(i, '-', gobbleID, ': has area:', newarea, 'which is bigger than the max:', gobbleareas_upper, '\n')
        
        # modifies the tree heights, this is actually only used for validation purposes and isn't really necessary.
          treeHeights[c(gobbleID,treeID)]<-max(treeHeights[c(gobbleID,treeID)])
          treeHeights[treeID]<-NA
        # for modifying the tree areas on the fly: 
          rightsize_gobbleID<-which(match_treeID)[combarea_test][which.max(G)] # a pointer to the position of the gobbling tree
                                                              # in the rightsize subset space.
          treeAreas_rightsize[rightsize_gobbleID]<-gobbleareas+treeAreas_toosmall[i]
        
        # adds the globbled area to the gobbling tree area.
        # modifies the tree areas:
          treeAreas[c(gobbleID,treeID)]<-treeAreas_rightsize[rightsize_gobbleID]
          treeAreas[treeID]<-NA
        }  
      }
    }
    treeHeights<-treeHeights[!is.na(treeHeights)]
  } 
  # End of gobble
  
  Crowns[coordCrown]<-tree.ind
  Index[coordSeeds]<-1:Ntrees
  
  if(cantext=='on')
  {
    coordSeeds_len<-nrow(coordSeeds)
    can_pix<-round(cansize/res(imagery)[1])
    row_max<-dim(Gnew)[1]
    col_max<-dim(Gnew)[2]
    subcan_mean<-vector('numeric', length=coordSeeds_len)
    subcan_cov<-vector('numeric', length=coordSeeds_len)
    subcan_range<-vector('numeric', length=coordSeeds_len)
    for(ind in 1:Ntrees)
    {
      prog.seq<-seq(0,1, by=0.1)
      progress<-round(coordSeeds_len * prog.seq)  %in% ind
      if(any(progress))
        cat(paste((which(progress)-1)*10), '%   ', sep='')
      
      row_pix<-coordSeeds[ind,'row']+c(-can_pix, can_pix)
      col_pix<-coordSeeds[ind,'col']+c(-can_pix, can_pix)
      in_img_test<-!(any(row_pix<=0 | row_pix>=row_max) | any(col_pix<=0 | col_pix>=col_max))
      if(in_img_test)
      {
        subcan<-Gnew[row_pix, col_pix]
        subcan_mean[ind]<-mean(subcan)
        subcan_sd<-sd(subcan)
        subcan_cov[ind]<-subcan_sd/subcan_mean[ind]
        subcan_range<-range(subcan)
        subcan_range[ind]<-subcan_range[2]-subcan_range[1]
      }
      else
      {
        subcan_mean[ind]<-NA
        subcan_cov[ind]<-NA
        subcan_range[ind]<-NA
      }  
    }
  }
  
  # converts back to a spatial grid:
  Cb <- imagery
  Mb <- imagery
  Cb[] <- as.numeric(Index[1:dim(Index)[1], dim(Index)[2]:1], 
                     byrow = TRUE) # given the correct orientation again.
  Cb[Cb==0]<-NA # Excludes non maxima pixels.
  
  m2 <- methods::as(Cb, "SpatialGridDataFrame")
  m3 <- raster::raster(m2, layer = 1)
  maxima.xy <- raster::rasterToPoints(m3, fun = , dissolve = TRUE)
  #IT <- sp::SpatialPointsDataFrame(m3.shp, data = data.frame(tree=m3.shp[,'layer']), match.ID = F)
  
  Cb <- imagery
  Mb <- imagery
  Cb[] <- as.numeric(Crowns[1:dim(Crowns)[1], dim(Crowns)[2]:1], 
                     byrow = TRUE)
  Mb[] <- as.numeric(Max[1:dim(Max)[1], dim(Max)[2]:1], 
                     byrow = TRUE)
  m2 <- methods::as(Cb, "SpatialGridDataFrame")
  m3 <- raster::raster(m2, layer = 1)
  
  m3.shp<-gdal_polygonizeR(m3, outshape=NULL, gdalformat = 'ESRI Shapefile',
                           pypath="/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_polygonize.py", readpoly=TRUE, quiet=TRUE)
  #m3.shp <- raster::rasterToPolygons(m3, fun = , dissolve = TRUE)
  
  # Extracts the mean and max tree heights:      
  names(m3.shp@data) <- "tree"
  m3.shp<-m3.shp[m3.shp$tree!=0,]
  treepos.ind<-match(1:Ntrees, m3.shp$tree)
  m3.shp[treepos.ind, "CH_mean"]<-treeHeights_mean
  m3.shp[treepos.ind, "CH_max"]<-treeHeights
  # The proportion of the total edge violations caused by each boundary condition:
  m3.shp[treepos.ind, "hbound"]<-boundcount[,1]
  m3.shp[treepos.ind, "sobelbound"]<-boundcount[,2]
  m3.shp[treepos.ind, "crownbound"]<-boundcount[,3]
  m3.shp[treepos.ind, "allombound"]<-boundcount[,4]
  m3.shp[treepos.ind, "tallerbound"]<-boundcount[,5]
  if(cantext=='on'){
    m3.shp[treepos.ind, "sobel"]<-sobel_mean
    m3.shp[treepos.ind, "can_mean"]<-subcan_mean
    m3.shp[treepos.ind, "can_cov"]<-subcan_cov
    m3.shp[treepos.ind, "can_range"]<-subcan_range
  }
  HCbuf <- rgeos::gBuffer(m3.shp, width = -res(imagery)[1]/2, byid = T) 
  ITCcv <- rgeos::gConvexHull(HCbuf, byid = T)
  ITCcvSD <- sp::SpatialPolygonsDataFrame(ITCcv, data = HCbuf@data, match.ID = F)
  #ITCcvSD <- sp::SpatialPolygonsDataFrame(HCbuf, data = HCbuf@data, match.ID = F)
  # Buffering seems silly... let's just delete the holes.
  
  #par(mfrow=c(1,2), mar=c(0,0,0,0))
  #plot(HCbuf)
  #plot(ITCcv)
  #plot(ITCcvSD)
  
  ITCcvSD$CA_m2 <- unlist(lapply(ITCcvSD@polygons, function(x) methods::slot(x, "area")))
  ITCcvSD$CR_m <- sqrt(ITCcvSD$CA_m2/pi)
  
  #par(mfrow=c(1,1), mar=c(4,4,2,2))
  #plot(y=ITCcvSD$CR_m, x=ITCcvSD$CH_max)
  
  # The xy location of the maxima:
  match.id<-match(ITCcvSD$tree, maxima.xy[,'layer'])
  ITCcvSD$maxx<-maxima.xy[match.id, 'x']
  ITCcvSD$maxy<-maxima.xy[match.id, 'y']
  
  #ITCcvSD <- ITCcvSD[ITCcvSD$CA_m2 > 1, ] # excludes the trees smaller than 1m^2 
  #return<-ITCcvSD
  
  #par(mfrow=c(1,1));plot(imagery);plot(ITCcvSD, add=TRUE, border='black')
  
  return(ITCcvSD)
}  