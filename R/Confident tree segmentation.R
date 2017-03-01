rm(list=ls())

# Header -----

library(maptools)
library(rgdal)
library(raster)
library(ggplot2)
library(gridExtra)
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
source("R/dtm_model_predict.R")
source("R/ggmultiplot.R")
source("R/DEMderiv.R")
source('R/min_extent.R')
source('R/poly_rast_vals.R')
source('R/crown_overlap.R')
source("R/rq_lut.R") # The function to run the quantile regression
source("R/Image preprocessing tools.R") # The function to run the quantile regression

LOGMODEL<-FALSE

# Load the UAV and LiDAR images ----
uav_dsm<-raster("data/uav_dsm_matched_cropped.tif")
uav_dtm<-raster("data/uav_dtm_matched_cropped.tif")
lid_dsm<-raster("data/lid_dsm_matched_cropped.tif")
lid_dtm<-raster("data/lid_dtm_matched_cropped.tif")
lid_chm<-raster("data/lid_chm_matched_cropped.tif")

lid_dsm<-lid_chm+lid_dtm # first a pit free LiDAR DSM is made from the CHM + DTM

# Measure the differences between the LiDAR and UAV values: -----

hist(values(lid_dsm))
hist(values(uav_dsm))

median(values(lid_dsm), na.rm=T)
median(values(uav_dsm), na.rm=T)

# ... there is approximately a 20 m difference between the elevations of the LiDAR and UAV
# DSMs. However by checking ground values in QGIS the difference is closer to 18 m. 


# Aligning the UAV and LiDAR models ----

# The thing to do is to see if you end up with negative or positive values on the roads etc once you subtract the 
# LiDAR CHM from the UAV CHM.

# Subtracts the difference between the two rasters from the UAV DSM:
uav_dsm<-uav_dsm-17.6
uav_dtm<-uav_dtm-17.6
lid_dsm<-lid_dsm-0 # Adding zero causes these objects to be stored in the memory which greatly speeds up operations later on.
lid_dtm<-lid_dtm-0
lid_chm<-lid_chm-0

# remove very low values from UAV DSM
uav_dsm[uav_dsm<0]<-NA
uav_dtm[uav_dtm<0]<-NA

# Calculates the canopy height models for the UAV
uav_chm<-uav_dsm-uav_dtm # UAV DTM

# Loads in the manual segmentation data ----

# Reads in the manually segmented trees:
man_uav<-readShapeSpatial("/Users/Tom/Documents/Work/R projects/UAVforestR/data/Manual trees UTM.shp")
man_lid<-readShapeSpatial("/Users/Tom/Documents/Work/R projects/UAVforestR/data/Manual trees LiDAR UTM.shp")


# Shift the manual segmentation according to a manual assessment of location against the LiDAR imagery
# in order to coregister them:
#manual<-shift(manual, x=-1, y=-4)
man_uav$R<-sqrt(man_uav$Area/pi) # Calculate a pseudoradius
man_lid$R<-sqrt(man_lid$Area/pi) # Calculate a pseudoradius

par(mfrow=c(1,1))
plot(uav_chm)
plot(man_uav, add=TRUE, border='red')
plot(man_lid, add=TRUE, border='blue')

# Calculate the height of the trees ----
man_uav$CH_mean<-poly_rast_vals(man_uav, uav_chm, FUN=mean)
man_uav$CH_max<-poly_rast_vals(man_uav, uav_chm, FUN=max)
man_lid$CH_mean<-poly_rast_vals(man_lid, lid_chm, FUN=mean)
man_lid$CH_max<-poly_rast_vals(man_lid, lid_chm, FUN=max)

# The relationship between mean and max height
ggplot(man_uav@data, aes(x=CH_max, y=CH_mean)) + geom_point() + geom_abline(intercept=0, slope=1)
ggplot(man_lid@data, aes(x=CH_max, y=CH_mean)) + geom_point() + geom_abline(intercept=0, slope=1)

# The allometric relationship
g1_uav<-ggplot(man_uav@data, aes(x=R, y=CH_max)) + geom_point() + geom_smooth(method='lm') + ylim(0,50)
g1_lid<-ggplot(man_lid@data, aes(x=R, y=CH_max)) + geom_point() + geom_smooth(method='lm') + ylim(0,50)
grid.arrange(g1_uav + ggtitle('UAV'), g1_lid + ggtitle('LiDAR'), ncol=2) 

par(mfrow=c(2,1));hist(man_uav$R, xlim=c(0,12)); hist(man_lid$R, xlim=c(0,12))

# Assessing the allometric relationship for the manually segmented trees ----

# Calculate the 90th percentile of the manual crowns height predicted from radius:
lut_lid<-rq_lut(x=man_lid$CH_max, y=man_lid$R, log=LOGMODEL)
lut_uav<-rq_lut(x=man_uav$CH_max, y=man_uav$R, log=LOGMODEL)
lut_rtoh_lid<-rq_lut(x=man_lid$R, y=man_lid$CH_max, log=LOGMODEL)
lut_rtoh_uav<-rq_lut(x=man_uav$R, y=man_uav$CH_max, log=LOGMODEL)

# Assessing the quantiles:
#par(mfrow=c(1,1))
#y<-allom_lookup(x, lut_rtoh_lid, tau=50)
#points(y~x, type='l')

# The mean difference between the LiDAR and UAV heights
mean(man_lid$CH_max-man_uav$CH_max)

# Crop to a smaller area
e1<-extent(313700, 313900, 9746700, 9746900)
uav_chm_e1<-crop(uav_chm, e1)
plot(uav_chm_e1)

# Height against radius 
# Using LiDAR allometry

par(mfrow=c(1,1), mar=c(4,4,2,2))
x<-seq(0,12,0.1)
plot(x, y=allom_lookup(x, lut_rtoh_lid, tau=90, antilog=LOGMODEL), type='l', ylab='Height', xlab='Radius', lwd=3, # 90th %ile
     ylim=c(0,50), xlim=c(0,10))
points(x, y=allom_lookup(x, lut_rtoh_lid, tau=50, antilog=LOGMODEL), type='l', col=rgb(0,0,0,1), lwd=3, lty=2) # 50th %ile
points(x, y=allom_lookup(x, lut_rtoh_lid, tau=20, antilog=LOGMODEL), type='l', col=rgb(0,0,0,1), lwd=3, lty=3) # 20th %ile
points(man_lid$R, y=man_lid$CH_max, col=rgb(0,0,0,0.5), pch=16) # LiDAR allometry
points(man_uav$R, y=man_uav$CH_max, col=rgb(0.7,0,0.3,0.5), pch=16) # UAV allometry 

# Load in the automatically segmented trees ----

# Load in the automatically segmented trees using the best set of parameters
# See the 'Testing all segmentation parameters' script for how this selection
# was made. The final choice was made according to the parameters that 
# had a size ratio and overlap between the automatic and manually segmented trees
# close to 1, and as many trees as possible with a shapebound score >0.9.
#

# Load in the data:
auto_uav<-readOGR("Data/ITC trees_params/seed_0.5_crown_0.7_sobel_3.5.shp")
match_uav<-readOGR("Data/ITC trees_params/Matched/seed_0.5_crown_0.7_sobel_3.5.shp")

# Exclude trees with a low shapebound score:
conf_match_uav<-match_uav[!is.na(match_uav$shapebound),]
conf_match_uav<-conf_match_uav[conf_match_uav$shapebound>0.9,]

# 
#conf_man_lid<-man_lid[man_lid$id %in% conf_match_uav$id,]
#lut_rtoh_conf_match_uav<-rq_lut(x=conf_man_lid$CR_m, y=conf_match_uav$CH_max)

# Exclude trees with a low shapebound score:
auto_uav@data$shapebound<-rowSums(auto_uav@data[,c('hbound', 'sobelbound')])
conf_auto_uav<-auto_uav[!is.na(auto_uav$shapebound),]
conf_auto_uav<-conf_auto_uav[conf_auto_uav$shapebound>0.9,]
nrow(conf_auto_uav)

# Plots the remaining data
par(mfrow=c(1,1));plot(uav_chm_e1)
plot(conf_auto_uav, add=TRUE)
plot(conf_match_uav, add=TRUE, col=rgb(0.8, 0, 0, 0.5))


# Plot the amount of overlap between the manual and autotrees:
par(mfrow=c(1,2))
hist(conf_match_uav$overlap_au, main='Proportion overlap', xlab=''); median(conf_match_uav$overlap_au)
hist(conf_match_uav$size_ratio, breaks=10, main='Size ratio', xlab=''); median(conf_match_uav$size_ratio)

par(mfrow=c(1,1))
plot(x=conf_auto_uav$CR_m, conf_auto_uav$CH_max, pch=16, col=rgb(0.7,0,0.3,0.5))
points(x=man_uav$R, man_uav$CH_max, pch=16, col=rgb(0,0,0,0.7))
points(x=man_lid$R, man_lid$CH_max, pch=16, col=rgb(0,0,1,0.7))
points(x, y=allom_lookup(x, lut_rtoh_lid, tau=50, antilog=LOGMODEL), type='l', col=rgb(0,0,0,1), lwd=3, lty=3) # 50th %ile

# Graphs showing the relationship between the overlap / size ratio and the other
# segmentation values:
ggplot(conf_match_uav@data, aes(x=CH_max, y=size_ratio)) + geom_point() + 
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, 3), se = FALSE)
ggplot(conf_match_uav@data, aes(x=CH_max, y=CR_m, size=size_ratio)) + geom_point()
ggplot(conf_match_uav@data, aes(x=CH_max, y=CR_m, size=overlap_au)) + geom_point()
ggplot(conf_match_uav@data, aes(x=CH_max, y=CR_m, size=size_ratio, color=shapebound)) + geom_point()

# Predict the size ratio:
#fm_ratio<-lm(size_ratio ~ splines::bs(CH_max, 3), data=conf_match_uav)
#summary(fm_ratio)
#ratio_pred<-predict(fm_ratio, newdata=data.frame(CH_max=conf_auto_uav$CH_max))
#par(mfrow=c(1,1));plot(ratio_pred~conf_auto_uav$CH_max)
#ratio_pred<-mean(conf_match_uav$size_ratio)
#CA_scaled<-(1/ratio_pred) * conf_auto_uav$CA_m2
#CR_m_scaled<-sqrt(CA_scaled/pi)
#conf_auto_uav$CH50<-allom_lookup(CR_m_scaled, lut_rtoh_lid, tau=68) 

# Estimate the tree height from the crown radius ----

# Convert the radius into height
conf_auto_uav$CH50<-allom_lookup(conf_auto_uav$CR_m, lut_rtoh_lid, tau=50, antilog=LOGMODEL)

# Currently this uses the 20th %ile... which doesn't seem right. It would be
# better to use the 50%ile I would think. I chose the 20th because it gives a good match
# to the LiDAR heights. 

# Checks the distribution of the supposed autotrees:
hist(conf_auto_uav$CR_m)
# Most of the trees have a crown radius less than 2 m but the allometric relationship
# is currently poorly defined at for trees of this size
# Having said that the error for these trees should be small 

# Extract the LiDAR chm values for the confidently segmented trees:
conf_auto_uav$lid_CH_max<-img_xy_extract(x=conf_auto_uav$maxx, y=conf_auto_uav$maxy, img=lid_chm)
# !!!! THIS SEEMS WELL DODGY - YOU COULD PERHAPS USE MEAN RATHER THAN MATCH THE HEIGHTS

par(mfrow=c(2,2), mar=c(4,4,2,2))
# Plots the UAV height against the LiDAR values:
plot(conf_auto_uav$CH_max~conf_auto_uav$lid_CH_max, pch=16, col=rgb(0,0,0,0.5), 
     xlab='LiDAR Height (m)', ylab='UAV Height (m)')
abline(0,1, col='blue')
# Plots the height estimates against the LiDAR values:
plot(conf_auto_uav$CH50~conf_auto_uav$lid_CH_max, pch=16, col=rgb(0,0,0,0.5), 
     xlab='LiDAR Height (m)', ylab='Height estimated from segmented crowns (m)', ylim=c(0,40))
abline(0,1, col='blue')
# Plots the height estimates again against the UAV heights:
plot(conf_auto_uav$CH50~conf_auto_uav$CH_max, pch=16, col=rgb(0,0,0,0.5), 
     xlab='UAV Height (m)', ylab='Height estimated from segmented crowns (m)')
points(conf_auto_uav$lid_CH_max~conf_auto_uav$CH_max, pch=16, col=rgb(0.7,0,0.3,0.5))
abline(0,1, col='blue')
# Plots the UAV and estimated heights against the segmented crown radii:
plot(conf_auto_uav$CH_max~conf_auto_uav$CR_m, pch=16, col=rgb(0,0,0,0.5), 
     xlab='Crown radius (m)', ylab='UAV Height (m)')
points(conf_auto_uav$CH50~conf_auto_uav$CR_m, pch=16, col=rgb(0.7,0,0.3,0.5))
points(conf_match_uav$CH_max~conf_match_uav$CR_m, pch=16, col=rgb(0,0,1,0.7))


# The median difference between the UAV heights and the 50%ile heights
par(mfrow=c(3,1), mar=c(4,4,2,2))
hist(conf_auto_uav$lid_CH_max-conf_auto_uav$CH_max, breaks=30, xlim=c(-10,20), main='Lidar - UAV')
median(conf_auto_uav$lid_CH_max-conf_auto_uav$CH_max, na.rm=TRUE)
hist(conf_auto_uav$CH50-conf_auto_uav$CH_max, breaks=30, xlim=c(-10,20), main='Segmented estimates - UAV')
median(conf_auto_uav$CH50-conf_auto_uav$CH_max)
# Pretty good I reckon! 
hist(conf_auto_uav$lid_CH_max-conf_auto_uav$CH50, breaks=30, xlim=c(-10,20), main='Lidar -Segmented estimates')
median(conf_auto_uav$lid_CH_max-conf_auto_uav$CH50, na.rm=TRUE)


# Writes the confidently segmented auto trees:
writeSpatialShape(conf_auto_uav, "Data/ITC trees_params/Confident_seed_0.5_crown_0.7_sobel_3.5.shp")

# Estimate the new DTM from the tree height estimates ----

# Produce a set of points for the estimated DTM at the tree bases according to CH50
# Estimate the new DTM surface
est_dtm<-dtm_predict(x=conf_auto_uav$maxx, y=conf_auto_uav$maxy, z=conf_auto_uav$CH50,
                  img=uav_dsm,
                  predict_grid_by=25,
                  k=1000)

col.pal<-colorRampPalette(c("white", "goldenrod1", "forestgreen", "limegreen"))( 10 )
col.breaks<-seq(20, 60, length=length(col.pal)+1)

par(mfrow=c(1,3))
plot(lid_dtm, col=col.pal, breaks=col.breaks, colNA='black', main='LiDAR')
plot(uav_dtm, col=col.pal, breaks=col.breaks, colNA='black', main='SFM raw')
plot(est_dtm, col=col.pal, breaks=col.breaks, colNA='black', main='SFM corrected')

par(mfrow=c(1,1))
plot(est_dtm, col=col.pal, breaks=col.breaks, colNA='black')
plot(conf_auto_uav, add=TRUE, border='red')

# Estimate the new CHM surface by subtracting the estimated DTM from the DSM.

est_chm<-uav_dsm-est_dtm
est_chm[est_chm<0]<-0

col.pal<-colorRampPalette(c("white", "goldenrod1", "forestgreen", "limegreen"))( 10 )
col.breaks<-seq(-1, 60, length=length(col.pal)+1)

par(mfrow=c(1,3))
plot(lid_chm, col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_chm, col=col.pal, breaks=col.breaks, colNA='black')
plot(est_chm, col=col.pal, breaks=col.breaks, colNA='black')

# Compare the LiDAR and estimated models ----

lid_chm_grid<-grid_fun(lid_chm, grid_by=50, fun=function(x) mean(x, na.rm=TRUE)) # for the LiDAR CHM
uav_chm_grid<-grid_fun(uav_chm, grid_by=50, fun=function(x) mean(x, na.rm=TRUE)) # for the LiDAR CHM
est_chm_grid<-grid_fun(est_chm, grid_by=50, fun=function(x) mean(x, na.rm=TRUE)) # for the LiDAR CHM

est_diff<-est_chm_grid-uav_chm_grid
hist(values(est_diff))

est_diff[est_diff<0]<-0
est_chm_grid2<-uav_chm_grid+est_diff

liduav_diff<-lid_chm_grid-uav_chm_grid
lidest_diff<-lid_chm_grid-est_chm_grid2

rmse<-function(x, y)  sqrt(mean((y-x)^2, na.rm=TRUE)) 

rmse(values(lid_chm_grid), values(uav_chm_grid))
rmse(values(lid_chm_grid), values(est_chm_grid))

col.pal<-colorRampPalette(c('red', "white", 'blue'))( 7 )
col.breaks<-seq(-14, 14, length=length(col.pal)+1)
par(mfrow=c(1,3), mar=c(2,2,2,2))
plot(liduav_diff, col=col.pal,  colNA="black", breaks=col.breaks, main="CHM difference")
plot(est_diff, col=col.pal, colNA="black", breaks=col.breaks, main="Correction amount")
plot(lidest_diff, col=col.pal, colNA="black", breaks=col.breaks, main="Corrected CHM difference")
#plot(conf_auto_uav, add=TRUE, border='black')


