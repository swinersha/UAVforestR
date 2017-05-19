
# There are 2 major problems with this approach:
# 1. We are massively undersegmenting big tree crowns because they have very strong edges and 
#    the crowns undulate more and so the minimum crown height increase rule seems to prevent 

# Change z_threshold and z_rise according to a measure of the position the maxima in the canopy
# Weight the mean by the strength of the edge.

require(raster)
require(dplyr)
require(tidyr)
require(ggplot2)

imagery<-raster("data/uav_dsm_matched.tif")
imagery<-crop(imagery, extent(313700, 313900, 9746600, 9746800))
plot(imagery)

distmax<-15
nlines<-12
lm.searchwin=3
profile_plot<-'off'
z_threshold<-3
z_rise<-4

source("R/net_make.R")

rtoh_lut<-read.csv("rtoh_lut_gta.csv")

rtoh<-function(x, lut, tau) 
{
  a<-lut[tau,'aa']
  b<-lut[tau,'bb']
  h<-(exp(a) * x^b)
  return(h)
}

#z<-linedat$h_diff; z_pairdiff<-linedat$z_pairdiff
zscale_edge<-function(z, z_pairdiff, z_threshold, z_rise)
{
  logic_fill<-function(x)
  {if(any(x)) x[which(x)[1]:length(x)]<-TRUE; return(x)}
  
  max_test<-logic_fill(z>0) # Any pixels after a height greater than the local maximum are excluded.
  diff_test<-logic_fill(z_pairdiff>z_rise) # Any pixesls after an increase of more than 20 cm are excluded; The increase should probably be scaled by the distance.
  #z_pairdiff<-z_pairdiff[!(max_test | diff_test)]
  z_pairdiff<-z_pairdiff[!diff_test]
  if(length(z_pairdiff)>=1)
  {
    if(min(z_pairdiff)<z_threshold) # If there are any differences greater than the threshold. Otherwise an NA value is returned.
    {
      edge_pos<-which.min(z_pairdiff) # Finds the stongest edge
      return(data.frame(pos=edge_pos, edge_strength=min(z_pairdiff)))
    }
    else
      return(data.frame(pos=NA, edge_strength=NA))
  }
  else
    return(data.frame(pos=NA, edge_strength=NA))
}  



#imagery=uav_chm
# Extracts just the locations of the tree maxima:
#itcIMG_fast<-function (imagery = NULL, TRESHSeed, 
#                       TRESHCrown, htod, specT, lm.searchwin=NULL, blur=TRUE, gobble='off')
#{
  # !!! You need to make some appropriate checks of the arguments here !!!
  cat('Blurring image\n')
  # Blurs the chm:
  if(blur)
  {
    imagery <- raster::focal(imagery, w = matrix(1, 3, 3), 
                             fun = function(x) {
                               mean(x, na.rm = T)
                             })
  }
  cat('Extracting matrix\n')
  # Extracts the image data as a matrix:
  #Max <- matrix(dim(imagery)[2], dim(imagery)[1], data = imagery[, 
                                                                 ], byrow = FALSE)
  cat('Preparing data structures\n')
  # Flips the image the right way again:
  #Max <- Max[1:dim(imagery)[2], dim(imagery)[1]:1]
  
  Max<-as.matrix(imagery)
  Gnew <- Max     # Copies the max matrix.
  Max[, ] <- 0    # Sets max to 0s. This will be used later for storing ...
  Index <- Max    # Copies max again.
  Index[, ] <- 0  # Sets index to 0s.
  Gnew[is.na(Gnew)] <- 0 # Sets any nas to 0s.
  #Gnew[Gnew < specT] <- 0 # Any values beneath the minimum height are set to 0s.
  index = 1       # Initiates the index
  II <- which(Gnew != 0, arr.ind = T) # Extracts the locations of pixels which are bigger than the min tree height.
  dim(II)
  
  #fm.win<-win.lmfun(win.min, win.max, z.min, z.max) # generates the relationship between height and search window size.
  # Extracts only the pixels that are sufficiently far from the image edge for the search window. In each direction.
  # The search window is selected according to the height of the tree:
  z<-Gnew[II[,1]+nrow(Gnew)*(II[,2]-1)] # extracts the tree heights from the matrix
  #z.pix<-z/res(imagery)[1] # converts heights in m to pixels
  #WinSize<-allom.dist(z.pix) # Finds the window size for each tree pixel in m
  #WinSize<-allom.dist(z, scale=res(imagery)[1], lut, tau=10)
  #pix2edge<-edge.dist(z, scale=res(imagery)[1], lut, tau=90) # Halfs the window size for subsetting by pixels far enough from the image edge
  
  # Extracts only the pixels far enough from the image edge.
  
  distmax_scale<-distmax/res(imagery)[1]
  II.ind<-II[, 1] > distmax_scale & 
    II[, 1] < (nrow(Gnew) - distmax_scale) &
    II[,2] > distmax_scale &
    II[,2] < (ncol(Gnew) - distmax_scale)
  II<-II[II.ind,]
  #WinSize<-WinSize[II.ind]
  z<-z[II.ind]
  dim(II)
  
  cat('Ordering pixels\n')
  # reorder from greatest to least according to z value (i.e. big trees first):
  z.order<-order(z, decreasing = TRUE)
  II<-II[z.order,]
  #WinSize<-WinSize[z.order]
  
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
  treeHeights<-vector('numeric', length=length.total)
  boundcount<-matrix(0, nrow=length.total, ncol=4, dimnames=list(NULL, c('crown2crown', 'allombound', 'heightbound', 'mheightbound')))
  tree.ind<-vector(length.total, mode='numeric')
  
  net<-net_make(nlines=nlines, distmax=distmax, im_res=res(imagery)) # makes the net
  
  net_xy<-net[,c('x', 'y')]
  colnames(net_xy)<-c('row', 'col')
  new_xy<-as.matrix(net_xy) # produces a fast indexing matrix for the net.
  
  net$h50<-rtoh(net$dist, rtoh_lut, 50) # calculates tree heights at the 50% for the crown diameter
  diff_ref60 <- -(net$h50 - net$h50 *0.6)
  diff_ref70 <- -(net$h50 - net$h50 *0.7)
  diff_ref80 <- -(net$h50 - net$h50 *0.8)
  diff_ref90 <- -(net$h50 - net$h50 *0.9)

  # Summary functions to be used in the loop below:
  CoV<-function(x) sd(x, na.rm=T) / mean(x, na.rm=T)
  NaN_to_NA<-function(x) ifelse(!is.nan(x), x, NA)
  
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
  
  # -----------------------------------------
  # Creates a matrix to store the summary of the distance estimation:
  # -----------------------------------------
  rad<-cov<-vector('numeric', length=Ntrees)
  nNA<-vector('integer', length=Ntrees)
  radius_est<-data.frame(rad, cov, nNA)
  TREE<-1
  plot(imagery)
  
  
  for (TREE in 1:Ntrees)
  {
  y<-coordSeeds[TREE,'row']
  x<-coordSeeds[TREE,'col']
  
  net_mat<-cbind(row=net_xy[,'row']+y, col=net_xy[,'col']+x)
  net_z<-Gnew[net_mat] # Extracts the height values.
  
  h_diff<-net_z-Gnew[y,x] # The difference in height at points on the line from the local maximum
  point_prof<-data.frame(line=net$line, x=net_mat[,'col'], y=net_mat[,'row'], dist=net$dist, pairdist=net$pairdist, z = net_z, h_diff = h_diff, diff_ref60 = diff_ref60, diff_ref70 = diff_ref70, diff_ref80 = diff_ref80, diff_ref90 = diff_ref90) # Makes a dataframe.
  point_prof$z_pairdiff<-unlist((lapply(1:nlines, function(LINE) c(0, diff(net_z[net$line==LINE]))))) # calculates the pairwise differences.
  
  can_pos<-mean(sapply(1:nlines, function(LINE) mean(point_prof$z_pairdiff[point_prof$line==LINE])))
  can_pos<--can_pos*distmax / res(imagery)[1] # The average z change across the net.
  
  adj_z_rise<-can_pos*0.1
  adj_z_rise[adj_z_rise<0]<-0
  
  if(profile_plot=='on')
  {
    g<-ggplot(point_prof, aes(x=dist, y=h_diff, col=as.factor(line))) + geom_point() + geom_line() + 
      geom_line(aes(dist, diff_ref60), color = 'black') +
      geom_line(aes(dist, diff_ref70), color = 'black') +
      geom_line(aes(dist, diff_ref80), color = 'black') +
      geom_line(aes(dist, diff_ref90), color = 'black')
    
    print(g)
    Sys.sleep(3)
  }

  # Test whether the pixels are taller than the maximum - if they are exclude them
  # Test whether there is a positive pixel height difference > 25 cm exclude any points after that
  # Find the edge as the greatest negative rate of change, on the condition that the rate must be greater
  # than some threshold (this could be scaled according to distance; but that could have its own problems)

  # Plots the pairwise differences in heights between the pixels:
  #ggplot(point_prof[point_prof$line==1,], aes(dist, z_pairdiff, color=as.factor(line))) +geom_line()
  #ggplot(point_prof, aes(dist, z_pairdiff, color=as.factor(line))) +geom_line()
   
  pp_edge<-vector('numeric', nlines)
  x_edge<-vector('numeric', nlines)
  y_edge<-vector('numeric', nlines)
  dist_edge<-vector('numeric', nlines)
  edge_strength<-vector('numeric', nlines)
  i<-4
  for(i in 1:nlines)
  {
    linedat<-point_prof[point_prof$line==i, ]
    edge_index<-with(linedat, zscale_edge(h_diff, z_pairdiff, z_threshold = z_threshold, z_rise = adj_z_rise))
    pp_edge[i]<-edge_index$pos
    if(!is.na(edge_index$pos))
    {
      x_edge[i]<-linedat$x[edge_index$pos]
      y_edge[i]<-linedat$y[edge_index$pos]
      dist_edge[i]<-linedat$dist[edge_index$pos]
      edge_strength[i]<-edge_index$edge_strength
    }
    else
    {
      x_edge[i]<-NA
      y_edge[i]<-NA
      dist_edge[i]<-NA
      edge_strength[i]<-NA
    } 
  }
  
  dist_edge
  
  radius_est[TREE, 1]<-mean(dist_edge, na.rm=TRUE)
  radius_est[TREE, 2]<-CoV(dist_edge)
  radius_est[TREE, 3]<-sum(is.na(dist_edge))

  # plot(imagery)
  # Create spatial polygons for each tree:
  if(nlines-radius_est[TREE,3]>5 & radius_est[TREE,2]<0.5)
  {
    tree.ind[count]<-TREE; count<-count+1
    angles <- seq(0,2*pi,length.out=360)
    point_circle<-data.frame(x=extent(imagery)@xmin + x *res(imagery)[1] + radius_est[TREE, 1]*cos(angles), y=extent(imagery)@ymax - y *res(imagery)[1] + radius_est[TREE, 1]*sin(angles))
    pp_xy[,'y']<-extent(imagery)@ymax - pp_xy$y *res(imagery)[1] 
    pp_xy[,'x']<-extent(imagery)@xmin + pp_xy$x *res(imagery)[1] 
    point_circle <- SpatialPolygons(list(Polygons(list(Polygon(point_circle)), ID=1)), proj4string=crs(imagery))
    point_circle_df <- SpatialPolygonsDataFrame(point_circle, data=data.frame(ID=1))
    plot(point_circle_df, add=TRUE)
    
    # Plots the rough polygons:
    #pp_xy<-data.frame(x=x_edge, y=y_edge)
    #pp_xy<-pp_xy[complete.cases(pp_xy),]
    #pp_xy<-pp_xy[c(1:nrow(pp_xy), 1), ]
    #pp_xy[,'y']<-extent(imagery)@ymax - pp_xy$y *res(imagery)[1] 
    #pp_xy[,'x']<-extent(imagery)@xmin + pp_xy$x *res(imagery)[1] 
    #sp_poly <- SpatialPolygons(list(Polygons(list(Polygon(pp_xy)), ID=1)), proj4string=crs(imagery))
    #sp_poly_df <- SpatialPolygonsDataFrame(sp_poly, data=data.frame(ID=1))
    #plot(sp_poly_df, add=TRUE)
  }  
}
 
  
  # -----------------------
  # Summary of the quality of the data.
   
  # Assesses the Usefulness of these data
  # In terms of the coefficient of variation at each height prop used:
  hprop_cov<-tidyr::gather(radius_est, key='hprop', value='cov', cov60:cov90, -c(dist60:dist90, nNA60:nNA90)) %>%
    select(hprop, cov)
  ggplot(hprop_cov, aes(cov, fill=hprop)) +geom_density(alpha = 0.3)
  # The variation in the estimates of crown radius decrease as the height proportion used
  # as the threshold for segmentation is decreased.
  
  select(radius_est, rad60:rad90) %>%
    summarise_each(funs(sum(!is.na(.))))
  # The number of trees also decreases with the height segmentation threshold.
  
  hprop_nNA<-tidyr::gather(radius_est, key='hprop', value='nNA', nNA60:nNA90) %>%
    select(hprop, nNA)
  ggplot(hprop_nNA, aes(nNA, fill=hprop)) +geom_density(alpha = 0.3)
  # A major reason for the low coefficients of variation when lower height thresholds are used
  # seems to be that there are far fewer measurements. In fact the majority come from only 1 measurement
  # So this doesn't seem very good.
  
  hprop_mean<-tidyr::gather(radius_est, key='hprop', value='mean', rad60:rad90) %>%
    select(hprop, mean)
  ggplot(hprop_mean, aes(mean, fill=hprop)) +geom_density(alpha = 0.3)
  # The radius estimates seem to be fairly similar irrespective of the threshold chosen. 
  # As such I would go with either 70% or 80% for the threshold. 
  # !!! These should be tested !!!
  # !!! We can than add an additional threshold to these to exclude trees with large COV values !!!
  dim(filter(hprop_cov, hprop=='cov80', cov<0.5))
  
  
  
  
  
  
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
  
  plot(m3.shp)
  #m3.shp <- raster::rasterToPolygons(m3, fun = , dissolve = TRUE)
  
  # Extracts the mean and max tree heights:      
  names(m3.shp@data) <- "tree"
  m3.shp<-m3.shp[m3.shp$tree!=0,]
  
  treepos.ind<-match(1:Ntrees, m3.shp$tree)
  m3.shp[treepos.ind, "CH_mean"]<-treeHeights_mean
  m3.shp[treepos.ind, "CH_max"]<-treeHeights
  # The proportion of the total edge violations caused by each boundary condition:
  m3.shp[treepos.ind, "crown2crown"]<-boundcount[,1]
  m3.shp[treepos.ind, "allombound"]<-boundcount[,2]
  m3.shp[treepos.ind, "heightbound"]<-boundcount[,3]
  m3.shp[treepos.ind, "mheightbound"]<-boundcount[,4]
  
  HCbuf <- rgeos::gBuffer(m3.shp, width = -res(imagery)[1]/2, byid = T) 
  ITCcv <- rgeos::gConvexHull(m3.shp, byid = T)
  plot(ITCcv)
  ITCcvSD <- sp::SpatialPolygonsDataFrame(ITCcv, data = HCbuf@data, match.ID = F)
  
  
  
  