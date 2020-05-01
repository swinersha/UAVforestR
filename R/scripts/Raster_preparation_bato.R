

rm(list=ls())

library(maptools)
library(rgdal)
library(raster)
library(ggplot2)
library(splancs)
library(RColorBrewer)
library(rgeos)
library(caret)

# Loads project functions:
#sapply(list.files(pattern="[.]R$", path="R/", full.names=TRUE), source)
source('R/min_extent.R')

lid_chm<-raster("data/raster/lid_chm_bato_50cm.tif")
uav_dsm<-raster("data/raster/Products/DSM_photoscan_wo_markers7and8.tif")
uav_dtm<-raster("data/raster/Products/DTM_photoscan_wo_markers7and8.tif")
uav_chm<-uav_dsm-uav_dtm

#------------------------------------------------

# Crop the lidar rasters to the extent of the UAV DSM
lid_chm<-crop(lid_chm, min_extent(uav_chm, lid_chm))
lid_dtm<-crop(lid_dtm, min_extent(uav_chm, lid_chm))
uav_chm<-crop(uav_chm, min_extent(uav_chm, lid_chm))

# Resample the UAV data so that they are on the same grid as the LiDAR data:
uav_chm_rs<-resample(uav_chm, lid_chm)
uav_dsm_rs<-resample(uav_dsm, lid_chm)


par(mfrow=c(1,3))
plot(lid_chm)
plot(uav_chm)
plot(uav_chm_rs)

writeRaster(uav_chm_rs, "data/raster/Products/CHM_photoscan_wo_markers7and8.tif")


# At this point I georeferenced the uav chm to the lidar using QGIS. The ground control points
# are saved in the products folder
# The resulting raster then gets read back and any values less than zero are
# converted to NA. This is specifcally because QGIS creates a load of really negative values.
# It is then realigned with the LiDAR CHM
uav_chm<-raster("data/raster/Products/CHM_photoscan_wo_markers7and8_georef.tif")
range(values(uav_chm), na.rm = TRUE)
uav_chm[uav_chm<0]<-NA
uav_chm[uav_chm>70]<-NA
uav_chm<-resample(uav_chm, lid_chm)

# Load the DSM too so that tpi can be calculated:
uav_dsm<-raster("data/raster/Products/DSM_photoscan_wo_markers7and8_georef.tif")
plot(uav_dsm)
range(values(uav_dsm), na.rm = TRUE)
uav_dsm<-resample(uav_dsm, lid_chm)


# Check the models against each other:
col.pal<-list(color = colorRampPalette(brewer.pal(9,"GnBu"))(10))$color
col.breaks<-seq(0, 50, length=length(col.pal)+1)
par(mfrow=c(1,2))
plot(uav_chm_rs, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
plot(uav_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)

# Check how well the models are aligned:
par(mfrow=c(1,2))
plot(lid_chm - uav_chm_rs) # Not well
plot(lid_chm - uav_chm) # Very well

aoi<-drawPoly(sp=TRUE, col='red', lwd=2)
plot(aoi, col = "blue", add = TRUE)
lid_chm<-mask(lid_chm, aoi)
uav_chm<-mask(uav_chm, aoi)
uav_dsm<-mask(uav_dsm, aoi)
uav_chm_rs<-mask(uav_chm_rs, aoi)
uav_dsm_rs<-mask(uav_dsm_rs, aoi)

lid_chm <- trim(lid_chm)
uav_chm <- trim(uav_chm)
uav_dsm<- trim(uav_dsm)
uav_chm_rs<- trim(uav_chm_rs)
uav_dsm_rs<- trim(uav_dsm_rs)

par(mfrow=c(1,1))
hist(lid_chm)

par(mfrow=c(1,2))
col.pal<-list(color = colorRampPalette(brewer.pal(9,"Greens"))(13))$color
col.breaks<-seq(0, 60, length=length(col.pal))
# LiDAR chm
plot(lid_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
# UAV chm
plot(uav_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)


writeRaster(lid_chm, "data/raster/lid_chm_bato_cropped.tif", overwrite = TRUE)
writeRaster(uav_chm, "data/raster/uav_chm_bato_georef_cropped.tif", overwrite = TRUE)
writeRaster(uav_dsm, "data/raster/uav_dsm_bato_georef_cropped.tif", overwrite = TRUE)
writeRaster(uav_chm_rs, "data/raster/uav_chm_bato_cropped.tif", overwrite = TRUE)
writeRaster(uav_dsm_rs, "data/raster/uav_dsm_bato_cropped.tif", overwrite = TRUE)


# Processing the bato dtm:
# This can be run alone.
lid_chm<-raster("data/raster/lid_chm_bato_cropped.tif")
lid_dtm<-raster("data/raster/lid_dtm_bato_1m.tif")
lid_dtm_rs<-resample(lid_dtm, lid_chm)

lid_dtm_rs<-crop(lid_dtm_rs, lid_chm)
lid_dtm_rs<-mask(lid_dtm_rs, lid_chm)

par(mfrow = c(1,2))
plot(lid_chm)
plot(lid_dtm_rs)

writeRaster(lid_dtm_rs, "data/raster/lid_dtm_bato_cropped.tif", overwrite = TRUE)
