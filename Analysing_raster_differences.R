rm(list=ls())

library(maptools)
library(rgdal)
library(raster)
library(ggplot2)
library(splancs)
library(RColorBrewer)
library(palettetown)
library(rgeos)
library(caret)
library(dplyr)
library(tidyr)
library(lmtest)
library(car)

source("R/build_dtm.R")
source("R/ITC segment update fast_gobble.R")
source("R/gdal_polygonizeR.R")
# Load the UAV and LiDAR images
uav_dsm<-raster("data/uav_dsm_matched.tif")
uav_dtm<-raster("data/uav_dtm_matched.tif")
lid_dsm<-raster("data/lid_dsm_matched.tif")
lid_dtm<-raster("data/lid_dtm_matched.tif")
lid_chm<-raster("data/lid_chm_matched.tif")

lid_dsm_solid<-lid_chm+lid_dtm # first a pit free LiDAR DSM is made from the CHM + DTM

n<-10
pokepal<-("Bulbasaur" %>% ichooseyou(5))[c(3,5)]
pokedex()
#------------------------------------------------

# Check for differences between the LiDAR and UAV values:

hist(values(lid_dsm_solid))
hist(values(uav_dsm))

median(values(lid_dsm_solid), na.rm=T)
median(values(uav_dsm), na.rm=T)

# ... there is approximately a 20 m difference between the elevations of the LiDAR and UAV
# DSMs. However by checking ground values in QGIS the difference is closer to 18 m. 
#------------------------------------------------

# The thing to do is to see if you end up with negative or positive values on the roads etc once you subtract the 
# LiDAR CHM from the UAV CHM.

# Subtracts the difference between the two rasters from the UAV DSM:
uav_dsm<-uav_dsm-16
uav_dtm<-uav_dtm-16
lid_dsm<-lid_dsm-0 # Adding zero causes these objects to be stored in the memory which greatly speeds up operations later on.
lid_dtm<-lid_dtm-0
lid_chm<-lid_chm-0

# remove very low values from UAV DSM
uav_dsm[uav_dsm<0]<-NA
uav_dtm[uav_dtm<0]<-NA
#------------------------------------------------
# Comparing the DSMs

dsm_diff<-uav_dsm-lid_dsm_solid # the difference is calculated
hist(values(dsm_diff))
median(values(dsm_diff), na.rm=T) # the median difference is 1 m
sqrt(mean(values(dsm_diff)^2, na.rm=T)) # the RMSE is 6 m

par(mfrow=c(1,3), mar=c(2,2,0,0)); plot(lid_dsm_solid, colNA='black'); plot(uav_dsm, colNA='black'); plot(dsm_diff, colNA='black')

#------------------------------------------------
# Comparing the DTMs

# Removes very low values from the DTM
uav_dtm[uav_dtm<10]<-NA

# Calculates the difference between the LiDAR and UAV DTMs:
dtm_diff<-uav_dtm-lid_dtm

# Plots the differences:
(dtm_diff.mean<-mean(values(dtm_diff), na.rm=T))
(dtm_diff.sd<-sd(values(dtm_diff), na.rm=T))
# As a histogram:
par(mfrow=c(1,1), mar=c(4,4,2,2));hist(values(dtm_diff)); abline(v=dtm_diff.mean, col='blue')
abline(v=dtm_diff.mean+dtm_diff.sd*1.96, col='red', lty=3); abline(v=dtm_diff.mean-dtm_diff.sd*1.96, col='red', lty=3)
# As terrain maps:
par(mfrow=c(1,3), mar=c(2,2,0,0)); plot(lid_dtm, colNA='black'); plot(uav_dtm, colNA='black'); plot(dtm_diff, colNA='black')

#------------------------------------------------
# Comparing the CHMs

# Calculates the canopy height models for the UAV
uav_chm<-uav_dsm-uav_dtm # UAV DTM
# Removes very low values from the UAV CHM: these are anomolous; check the uav_chm prior to this step if you want to check.
# plot(uav_chm)
uav_chm[uav_chm<min(values(lid_chm), na.rm=TRUE)]<-NA

# calculates the difference between the uav and lidar CHMs:
chm_diff<-uav_chm-lid_chm

# Looking at the numeric differences
par(mfrow=c(1,1), mar=c(4,4,2,2)); hist(chm_diff); median(values(chm_diff), na.rm=TRUE); mean(values(chm_diff), na.rm=TRUE)

par(mfrow=c(1,3), mar=c(2,2,0,0)); plot(lid_chm, colNA='black'); plot(uav_chm, colNA='black'); plot(chm_diff, colNA='black')
#zoom(chm_diff)

#------------------------------------------------
# Comparing the differences in each layer (DSM, DTM, CHM)

par(mfrow=c(1,3), mar=c(2,2,0,0));plot(dsm_diff, colNA='black'); plot(dtm_diff, colNA='black'); plot(chm_diff, colNA='black')

#------------------------------------------------
# Mean values differences by gridding the rasters

uav_dsm_grid<-grid_mean(uav_dsm, grid_by=25)
uav_dtm_grid<-grid_mean(uav_dtm, grid_by=25)
uav_chm_grid<-grid_mean(uav_chm, grid_by=25)
lid_dsm_grid<-grid_mean(lid_dsm, grid_by=25)
lid_dtm_grid<-grid_mean(lid_dtm, grid_by=25)
lid_chm_grid<-grid_mean(lid_chm, grid_by=25)

chm_grid_diff<-uav_chm_grid-lid_chm_grid # Calculates the difference between the UAV and lidar CHMs when gridded
chm_grid_diff[chm_grid_diff>5]<-NA #
cov_chm<-chm_grid_diff/lid_chm_grid # The relative size of the error
cov_chm[cov_chm>1]<-NA

# Plots out the gridded differnces
col.pal<-colorRampPalette(c(pokepal[1], "white", pokepal[2]))( 10 )
col.breaks<-seq(-5, 15, length=length(col.pal)+1)
par(mfrow=c(1,1), mar=c(2,2,2,2))
par(mfrow=c(1,3), mar=c(2,2,2,2)); plot(uav_dsm_grid-lid_dsm_grid, col=col.pal,  colNA="black", breaks=col.breaks, main="DSM")
plot(uav_dtm_grid-lid_dtm_grid, col=col.pal, colNA="black", breaks=col.breaks, main="DTM")
plot(lid_chm_grid-uav_chm_grid, col=col.pal, colNA="black", breaks=col.breaks, main='CHM')

# Another couple of plots
par(mfrow=c(1,4), mar=c(2,2,0,0)); plot(lid_dtm); plot(uav_dtm); plot(chm_grid_diff); plot(chm_diff)
par(mfrow=c(1,4), mar=c(2,2,0,0)); plot(lid_chm); plot(uav_chm); plot(chm_grid_diff); plot(cov_chm)

# The average difference between the UAV and LiDAR CHMs
par(mfrow=c(1,1)); hist(values(chm_grid_diff))
mean(values(chm_grid_diff), na.rm=TRUE); sd(values(chm_grid_diff), na.rm=TRUE) 


#------------------------------------------------
# Taking a segmentation approach
#------------------------------------------------

rtoh_lut<-read.csv("rtoh_lut_gta.csv")

# Calculates the height 50%ile for the observed crown radius:
rtoh<-function(x, lut, tau) 
{
  a<-lut[tau,'aa']
  b<-lut[tau,'bb']
  h<-(exp(a) * x^b)
  return(h)
}

#-------/\--/\--------
# ------'----'--------
# Working with a cropped area

extent(uav_chm)
e1<-extent(313450, 314025, 9746500, 9747000)
#e2<-extent(313500, 314500, 9746000, 9746500)
#e3<-extent(313000, 314000, 9747000, 9747500)

lid_dsm_crop<-crop(lid_dsm, e1)
uav_dsm_crop<-crop(uav_dsm, e1)
lid_dtm_crop<-crop(lid_dtm, e1)
uav_dtm_crop<-crop(uav_dtm, e1)
uav_chm_crop<-crop(uav_chm, e1)
lid_chm_crop<-crop(lid_chm, e1)

CH_lim<-5
hbound_lim<-0.4
sobel_lim<--5 # This is needed because sometimes the sobel value comes out as -Inf

htor_lut_gta<-read.csv("htor_lut_gta.csv")
lut<-htor_lut_gta
htod<-htod_hrf
uav_trees<-itcIMG_fast(uav_chm_crop+5, THRESHSeed=0.7, THRESHCrown=0.7, htod=htod_hrf, specT=2, SOBELstr=2, lm.searchwin=NULL, blur=TRUE, gobble='off')
uav_trees$normsobel<-log(uav_trees$sobel/uav_trees$CH_max)
uav_trees$logsobel<-log(uav_trees$sobel)
# Exploratory plots of the sobel edge intensity by tree height:
hist(uav_trees$sobelbound)
hist(uav_trees$hbound)
hist(uav_trees$logsobel)
par(mfrow=c(1,2), mar=c(4,4,2,2)); plot(y=uav_trees$logsobel, x=uav_trees$CH_max); plot(y=uav_trees$normsobel, x=uav_trees$CH_max)
# Plotting the CHMs
#plot(lid_chm_crop, colNA='black')
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(uav_chm_crop, colNA='black')
#plot(uav_trees, add=TRUE)
conf.CD<-uav_trees[which((uav_trees$hbound+uav_trees$sobelbound)>0.4 & 
                           uav_trees$CH_max>CH_lim),]

#conf.CD<-uav_trees[which((uav_trees$sobelbound)>0.4 & 
#                           uav_trees$CH_max>CH_lim),]
col.pal<-colorRampPalette(c("white", pokepal[2]))( 10 )
#col.pal<-c(col.pal, 'blue')
#col.breaks<-c(seq(0, 1, length=length(col.pal)+1), 5)
#cut_cols<-cut(conf.CD$normsobel, col.breaks)
col.breaks<-seq(-2, 4, length=length(col.pal)+1)
cut_cols<-cut(conf.CD$logsobel, col.breaks)
#cut_cols<-cut(uav_trees$heightbound+uav_trees$heightbound, col.breaks)
cut_cols<-col.pal[cut_cols]
#plot(uav_trees, add=TRUE, border=cut_cols)
plot(conf.CD, add=TRUE, border=cut_cols)

#plot(1:length(col.pal), col=col.pal, pch=16, cex=10)
plot(uav_trees$logsobel~uav_trees$CH_max); abline(h=20, col='red')
par(mar=c(4,4,2,2)); plot(y=uav_trees$sobel, x=(uav_trees$hbound+uav_trees$sobelbound)); abline(h=0, col='red')

#------------------------------------------------

# Set the tree segmentation to use in fitting the GAM:
conf.trees<-uav_trees

conf.trees$CH50<-rtoh(conf.trees$CR_m, rtoh_lut, 50)
hist(conf.trees$CH50)

# Extracts the DSM values for the confidently segmented trees:
treepos_xy<-data.frame(x=conf.trees$maxx, y=conf.trees$maxy)
treepos<-SpatialPoints(treepos_xy, proj4string = crs(uav_dsm))
treepos_dsm<-raster::extract(uav_dsm, treepos)
treepos_uavchm<-raster::extract(uav_chm, treepos)
treepos_lidchm<-raster::extract(lid_chm, treepos)

segtree<-data.frame(x=conf.trees$maxx, y=conf.trees$maxy, h_uav=treepos_uavchm, h_lid=treepos_lidchm, CR_m=conf.trees$CR_m, h_CH50=conf.trees$CH50, sobel=conf.trees$logsobel, logsobel=conf.trees$logsobel, normsobel=conf.trees$normsobel, hbound=conf.trees$hbound+conf.trees$sobelbound, sobelbound=conf.trees$sobelbound)

ggplot(segtree, aes(x=h_uav, y=h_lid, color=hbound)) + geom_point() + geom_smooth()
ggplot(segtree, aes(x=h_uav, y=h_lid, color=logsobel)) + geom_point() + geom_smooth()
ggplot(segtree, aes(x=h_uav, y=h_lid, color=CR_m)) + geom_point() + geom_smooth()
ggplot(segtree[segtree$hbound>0,], aes(x=h_uav, y=h_lid, color=normsobel)) + geom_point() + geom_smooth()
# It looks like there are some features which might be useful in 
# separating true 0 CHM values from non-0 CHM values.
ggplot(segtree[segtree$hbound>0.4,], aes(x=h_uav, y=h_lid, color=normsobel, size=CR_m)) + geom_point() + geom_smooth()
ggplot(segtree, aes(x=h_uav, y=CR_m, col=hbound)) + geom_point() + geom_smooth()
ggplot(segtree, aes(x=h_lid, y=CR_m, col=hbound)) + geom_point() + geom_smooth()
ggplot(segtree, aes(x=h_lid, y=CR_m, col=normsobel)) + geom_point() + geom_smooth()
ggplot(segtree[segtree$hbound>0.4,], aes(x=h_uav, y=CR_m, col=hbound)) + geom_point() + geom_smooth()
ggplot(segtree[segtree$hbound>0.4,], aes(x=h_lid, y=CR_m, col=hbound)) + geom_point() + geom_smooth(method = 'lm', formula=y ~ poly(x, 3, raw=TRUE))
# The relationship between the sobel edge, height and radius:
ggplot(segtree[segtree$hbound>0.4,], aes(x=h_uav, y=logsobel)) + geom_point() + geom_smooth()
ggplot(segtree[segtree$hbound>0.4,], aes(x=CR_m, y=logsobel)) + geom_point() + geom_smooth()


segtree_sub<-segtree[segtree$hbound>hbound_lim,]
segtree_sub<-segtree[segtree$sobel>sobel_lim,]
segtree_sub<-segtree_sub[complete.cases(segtree_sub),]
fm_treeest0<-lm(h_lid~h_uav, data=segtree_sub)
fm_treeest1<-lm(h_lid~h_uav+logsobel+hbound, data=segtree_sub)
fm_treeest2<-lm(h_lid~h_uav+CR_m*logsobel+hbound, data=segtree_sub)
fm_treeest3<-lm(h_lid~h_uav*CR_m*logsobel+hbound, data=segtree_sub)
fm_treeest4<-lm(h_lid~h_uav*CR_m+logsobel, data=segtree_sub)
fm_treeest5<-lm(h_lid~h_uav+CR_m+logsobel+h_uav:CR_m+h_uav:logsobel, data=segtree_sub)
fm_treeest6<-lm(h_lid~h_uav+CR_m+logsobel+h_uav:CR_m+CR_m:logsobel, data=segtree_sub)
fm_treeest7<-lm(h_lid~h_uav*CR_m*logsobel, data=segtree_sub)

#fm_gam1<-mgcv::gam(h_lid~s(h_uav)+s(x, y, k=500), data=segtree_sub)
#fm_gam2<-mgcv::gam(h_lid~s(h_uav, k=500)+s(x, y, k=500), data=segtree_sub)
#fm_gam3<-mgcv::gam(h_lid~s(h_uav, k=500)+CR_m+sobel+h_uav:CR_m+CR_m:sobel+s(x, y, k=500), data=segtree_sub) # having smoothing on UAV height is silly and will overfit
#fm_gam4<-mgcv::gam(h_lid~s(h_uav)+s(CR_m)+s(sobel)+h_uav:CR_m+CR_m:sobel+s(x, y, k=500), data=segtree_sub)
#fm_gam5<-mgcv::gam(h_lid~s(h_uav)+s(CR_m)+s(sobel)+s(x, y, k=500), data=segtree_sub)
#fm_gam6<-mgcv::gam(h_lid~s(h_uav)+s(CR_m)+sobel+s(x, y, k=500), data=segtree_sub) # All the above  models have smoothing terms which will probably lead to overfitting
# ... ---
#fm_gam7<-mgcv::gam(h_lid~h_uav+CR_m+sobel+s(x, y, k=500), data=segtree_sub)
#fm_gam8<-mgcv::gam(h_lid~h_uav+CR_m*sobel+s(x, y, k=500), data=segtree_sub)
#fm_gam9<-mgcv::gam(h_lid~h_uav+CR_m*sobel+s(x, y, k=500), data=segtree_sub, weights=(hbound/mean(hbound))) # the weights don't seem much good.
#fm_gam10<-mgcv::gam(h_lid~h_uav+CR_m*sobel+hbound+s(x, y, k=500), data=segtree_sub) # Adding hbound doesn't improve on gm_gam8 without it

# Test model fit
#AIC(fm_treeest0, fm_treeest1, fm_treeest2, fm_treeest3, fm_treeest4, fm_treeest5, fm_treeest6, fm_treeest7,
#    fm_gam1, fm_gam2, fm_gam3, fm_gam4, fm_gam5, fm_gam6, fm_gam7, fm_gam8, fm_gam9, fm_gam10)
AIC(fm_treeest0, fm_treeest1, fm_treeest2, fm_treeest3, fm_treeest4, fm_treeest5, fm_treeest6, fm_treeest7)

model<-fm_treeest3
summary(model)
par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(model)

segtree_sub$pred<-fitted(model)

# plot model predictions against reality:
par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(x=segtree$h_uav, y=segtree$h_lid)
points(x=segtree_sub$h_uav, y=segtree_sub$h_lid, col='red')
plot(x=fitted(model), y=segtree_sub$h_lid, xlim=c(0, max(segtree_sub$pred, na.rm=TRUE)))

ggplot(segtree_sub, aes(x=pred, y=h_lid, color=logsobel, size=CR_m)) + geom_point() + geom_smooth()
# The smallest tree we have predicted the height of is 10 m tall; this is because 5m was added to the CHM and then
# trees smaller than 5 m were excluded (see above).

# Plots out the predictions against the orginal UAV CHM values:
par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(x=segtree_sub$h_uav, y=segtree_sub$pred, xlim=c(0, max(segtree_sub$pred))); abline(0,1, col='red')
hist(segtree_sub$h_uav-segtree_sub$pred)

# Extracts the DSM values for the confidently segmented trees:
treepos_xy<-data.frame(x=segtree_sub$x, y=segtree_sub$y)
treepos<-SpatialPoints(treepos_xy, proj4string = crs(uav_dsm))
treepos_dsm<-raster::extract(uav_dsm, treepos)
treepos_uavchm<-raster::extract(uav_chm, treepos)
treepos_lidchm<-raster::extract(lid_chm, treepos)
# Calculates the estimated ground positions:
treepos_dtm<-treepos_dsm-segtree_sub$pred
treepos_dtm<-cbind(treepos_xy, z=treepos_dtm)

predict_grid_by<-5 # The resolution of the predicted DTM
imagery<-uav_chm_crop # The image to use in the prediction.
k<-100

# Estimating the ground 
cat("\nFitting GAM\n")
est.dtm <- mgcv::gam(z ~ s(x, y, k=k), data=treepos_dtm) # fit the 2D spline.
summary(est.dtm)
# Create a raster with the predicted values
grid<-img_grid(imagery, grid_size=predict_grid_by) # produce a grid for prediction.
dtm_pred<-mgcv::predict.gam(est.dtm, newdata=grid) # predict
dtm_pred<-cbind(grid, dtm_pred)
sp::coordinates(dtm_pred)=~x+y # convert to gridded spatial object
sp::proj4string(dtm_pred)<-raster::crs(imagery)
sp::gridded(dtm_pred) <-TRUE
dtm_pred <- raster::raster(dtm_pred) # convert to raster
dtm_pred<-resample(dtm_pred, imagery)
dtm_pred<-mask(dtm_pred, imagery)

chm_pred<-uav_dsm-dtm_pred

col.pal<-colorRampPalette(c("white", "goldenrod1", "forestgreen", "limegreen"))( 10 )
col.breaks<-seq(0, 60, length=length(col.pal)+1)

par(mfrow=c(2,3), mar=c(2,2,2,2))
plot(lid_dtm_crop)
plot(uav_dtm_crop)
plot(dtm_pred)
plot(crop(lid_chm, chm_pred), col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_chm_crop, col=col.pal, breaks=col.breaks, colNA='black')
plot(chm_pred, col=col.pal, breaks=col.breaks, colNA='black')


#--------------------------------------------
# Trying it on a naive area - for cross-validation:
#--------------------------------------------

extent(uav_chm)
#e1<-extent(313500, 314000, 9746500, 9747000)
e2<-extent(313500, 314500, 9746000, 9746500)
#e2<-extent(313000, 314000, 9747000, 9747500)

lid_dsm_crop<-crop(lid_dtm, e2)
uav_dsm_crop<-crop(uav_dtm, e2)
lid_dtm_crop<-crop(lid_dtm, e2)
uav_dtm_crop<-crop(uav_dtm, e2)
uav_chm_crop<-crop(uav_chm, e2)
lid_chm_crop<-crop(lid_chm, e2)

par(mfrow=c(1,1))
plot(uav_chm_crop)

uav_trees2<-itcIMG_fast(uav_chm_crop+5, THRESHSeed=0.7, THRESHCrown=0.7, htod=htod_hrf, specT=2, SOBELstr = 2, lm.searchwin=NULL, blur=TRUE, gobble='off')
uav_trees2$normsobel<-log(uav_trees2$sobel/uav_trees2$CH_max)
uav_trees2$logsobel<-log(uav_trees2$sobel)

# Extracts only the confidently segmented trees using the same parameters as
# when fitting the model to the LiDAR data originally:
conf.trees<-uav_trees2[which((uav_trees2$hbound+uav_trees2$sobelbound)>hbound_lim & 
                           uav_trees2$CH_max>CH_lim),]

# Extracts the DSM values for the confidently segmented trees:
treepos_xy<-data.frame(x=conf.trees$maxx, y=conf.trees$maxy)
treepos<-SpatialPoints(treepos_xy, proj4string = crs(uav_dsm))
treepos_dsm<-raster::extract(uav_dsm, treepos)
treepos_uavchm<-raster::extract(uav_chm, treepos)
treepos_lidchm<-raster::extract(lid_chm, treepos)

# Creates the prediction dataframe from the segmented trees:
segtree<-data.frame(x=conf.trees$maxx, y=conf.trees$maxy, h_uav=treepos_uavchm, h_lid=treepos_lidchm, CR_m=conf.trees$CR_m, sobel=conf.trees$sobel, logsobel=conf.trees$logsobel, normsobel=conf.trees$normsobel, hbound=conf.trees$hbound+conf.trees$sobelbound)
segtree_sub<-segtree[segtree$hbound>hbound_lim,]
segtree_sub<-segtree_sub[complete.cases(segtree_sub),]

segtree_sub$pred<-predict(model, newdata=segtree_sub)

# Plots out the GAM predictions:
par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(x=segtree_sub$h_uav, y=segtree_sub$pred); abline(0,1, col='red')
hist(segtree_sub$h_uav-segtree_sub$pred)

# Extracts the DSM values for the confidently segmented trees:
treepos_xy<-data.frame(x=segtree_sub$x, y=segtree_sub$y)
treepos<-SpatialPoints(treepos_xy, proj4string = crs(uav_dsm))
treepos_dsm<-raster::extract(uav_dsm, treepos)
treepos_uavchm<-raster::extract(uav_chm, treepos)
treepos_lidchm<-raster::extract(lid_chm, treepos)
treepos_dtm<-treepos_dsm-segtree_sub$pred
treepos_dtm<-cbind(treepos_xy, z=treepos_dtm, w=segtree_sub$logsobel) # You can try playing arond with what you weight by ...
                                          # I have tried h_uav, hbound and sobel, they are all pretty similar 
                                          # but sobel looks best

predict_grid_by<-5 # The resolution of the predicted DTM
imagery<-uav_chm_crop # The image to use in the prediction.
k=100 # Seetng this too high makes the model too bumpy; It might be worth playing with it though

# Estimating the ground 
cat("\nFitting GAM\n")
est.dtm <- mgcv::gam(z ~ s(x, y, k=k), data=treepos_dtm) #, weights=w/mean(w)) # fit the 2D spline.
summary(est.dtm)
# Create a raster with the predicted values
grid<-img_grid(imagery, grid_size=predict_grid_by) # produce a grid for prediction.
dtm_pred<-mgcv::predict.gam(est.dtm, newdata=grid) # predict
dtm_pred<-cbind(grid, dtm_pred)
sp::coordinates(dtm_pred)=~x+y # convert to gridded spatial object
sp::proj4string(dtm_pred)<-raster::crs(imagery)
sp::gridded(dtm_pred) <-TRUE
dtm_pred <- raster::raster(dtm_pred) # convert to raster
dtm_pred<-resample(dtm_pred, imagery)
dtm_pred<-mask(dtm_pred, imagery)

# Create a raster with the standard errors:
dtm_se<-mgcv::predict.gam(est.dtm, newdata=grid, se.fit=TRUE)$se.fit*2
dtm_se<-cbind(grid, dtm_se)
sp::coordinates(dtm_se)=~x+y # convert to gridded spatial object
sp::proj4string(dtm_se)<-raster::crs(imagery)
sp::gridded(dtm_se) <-TRUE
dtm_se <- raster::raster(dtm_se) # convert to raster
dtm_se<-resample(dtm_se, imagery)
dtm_se<-mask(dtm_se, imagery)
plot(dtm_se)

chm_pred<-uav_dsm-dtm_pred
uav_trees2_conf<-uav_trees2[which((uav_trees2$hbound+uav_trees2$sobelbound)>0.4),]

col.pal<-colorRampPalette(c("white", "goldenrod1", "forestgreen", "limegreen"))( 10 )
par(mfrow=c(2,3), mar=c(2,2,2,2))
col.breaks<-seq(20, 50, length=length(col.pal)+1)
plot(lid_dtm_crop, col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_dtm_crop, col=col.pal, breaks=col.breaks, colNA='black')
plot(dtm_pred, col=col.pal, breaks=col.breaks, colNA='black')
col.breaks<-seq(0, 60, length=length(col.pal)+1)
plot(crop(lid_chm, chm_pred), col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_chm_crop, col=col.pal, breaks=col.breaks, colNA='black')
plot(chm_pred, col=col.pal, breaks=col.breaks, colNA='black')

par(mfrow=c(1,1), mar=c(2,2,2,2))
plot(chm_pred, col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_trees2_conf, add=TRUE, border='blue')

# It seems to work!

#--------------------
# Checking again with the gridded values:

lid_chm_crop_grid<-grid_mean(lid_chm_crop, grid_by=25) # for the LiDAR CHM
uav_chm_crop_grid<-grid_mean(uav_chm_crop, grid_by=25) # for the new predicted CHM
chm_pred_grid<-grid_mean(chm_pred, grid_by=25) # for the new predicted CHM


par(mfrow=c(1,2))
plot(lid_chm_crop)
plot(chm_pred)

col.breaks<-seq(00, 30, length=length(col.pal)+1)
par(mfrow=c(1,2))
plot(lid_chm_crop_grid, col=col.pal, breaks=col.breaks, colNA='black')
plot(chm_pred_grid, col=col.pal, breaks=col.breaks, colNA='black')

chm_grid_diff<-lid_chm_crop_grid-uav_chm_crop_grid # Calculates the difference between the UAV and lidar CHMs when gridded
chm_pred_grid_diff<-lid_chm_crop_grid-chm_pred_grid # Calculates the difference between the UAV and lidar CHMs when gridded

# Calculat the RMSE:
sqrt(mean((values(lid_chm_crop_grid)-values(uav_chm_crop_grid))^2, na.rm=TRUE)) # Photoscan DTM
sqrt(mean((values(lid_chm_crop_grid)-values(chm_pred_grid))^2, na.rm=TRUE)) # PEstimated DTM.

# Plots this:
col.pal<-colorRampPalette(c(pokepal[1], "white", pokepal[2]))( 10 )
col.breaks<-seq(-12, 12, length=length(col.pal)+1)
par(mfrow=c(1,3), mar=c(2,2,2,2))
plot(chm_grid_diff, col=col.pal, colNA="black", breaks=col.breaks, main='UAV CHM error')
plot(chm_pred_grid_diff, col=col.pal, colNA="black", breaks=col.breaks, main='CHM prediction error')
col.pal<-colorRampPalette(c("white", "goldenrod1", "forestgreen", "limegreen"))( 10 )
col.breaks<-seq(0, 60, length=length(col.pal)+1)
plot(chm_pred, col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_trees2_conf, add=TRUE, border='blue')










#--------------------------------------------
# Calculates the DTM positions from the trees:
#--------------------------------------------

#------------------------------------------------

uav_dtm_smooth<-build_dtm(imagery=dsm, grid_by=50, k=100, predict_grid_by=10)
uav_dtm_smooth<-resample(uav_dtm_smooth, lid_dtm, filename="data/UAV DTM 33 cm UTM smooth.tif")

par(mfrow=c(1,3)); plot(lid_dtm, colNA='black'); plot(uav_dtm, colNA='black'); plot(uav_dtm_smooth, colNA='black')

dtmpred_diff<-uav_dtm_pred-lid_dtm
dtm_diff<-uav_dtm-lid_dtm
dtmpred_diff[dtmpred_diff<=(-15)]<-NA
par(mfrow=c(1,4)); plot(lid_dtm); plot(uav_dtm_pred); plot(dtmpred_diff); plot(dsm); plot(extent(dsm), add=TRUE)
par(mfrow=c(1,3)); plot(lid_dtm); plot(uav_dtm_pred); plot(uav_dtm); plot(dtm_diff)

#par(mfrow=c(1,3)); plot(lid_dtm); plot(uav_dtm); vis.gam(m.dtm, plot.type='contour')
#plot(dtm_pred)
#par(mfrow=c(1,3)); plot(lid_dtm); plot(uav_dtm); plot(dtm_pred)


hist(treepos_chm-conf.trees$CH50); mean(treepos_chm-conf.trees$CH50, na.rm=TRUE)
treepos_dtm<-treepos_dsm-conf.trees$CH50
hist(conf.trees$CH_max-conf.trees$CH50); median(conf.trees$CH_max-conf.trees$CH50, na.rm=TRUE)

treepos_dtm<-cbind(treepos_xy, z=treepos_dtm)

k<-1000
predict_grid_by<-20
imagery<-uav_dsm

cat("\nFitting GAM\n")
est.dtm <- mgcv::gam(z ~ s(x, y, k=k), data=treepos_dtm) # fit the 2D spline.
summary(est.dtm)
# Create a raster with the predicted values
grid<-img_grid(imagery, grid_size=predict_grid_by) # produce a grid for prediction.
dtm_pred<-mgcv::predict.gam(est.dtm, newdata=grid) # predict
dtm_pred<-cbind(grid, dtm_pred)
sp::coordinates(dtm_pred)=~x+y # convert to gridded spatial object
sp::proj4string(dtm_pred)<-raster::crs(imagery)
sp::gridded(dtm_pred) <-TRUE
dtm_pred <- raster::raster(dtm_pred) # convert to raster
dtm_pred<-resample(dtm_pred, imagery)
dtm_pred<-mask(dtm_pred, imagery)


col.pal<-colorRampPalette(c("white", "goldenrod1", "forestgreen", "limegreen"))( 10 )
col.breaks<-seq(10, 60, length=length(col.pal)+1)

par(mfrow=c(1,3), mar=c(2,2,2,2))
plot(lid_dtm, colNA='black', col=col.pal, breaks=col.breaks)
plot(uav_dtm, colNA='black', col=col.pal, breaks=col.breaks)
plot(dtm_pred, colNA='black', col=col.pal, breaks=col.breaks)



