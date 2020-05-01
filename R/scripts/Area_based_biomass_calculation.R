library(mgcv)
library(UAVforestR)
library(gridExtra)
library(sp)
library(rgeos)
library(rgdal)
library(RColorBrewer)
library(spatialEco)

htod<-function(x, tau) htod_lookup(x, lut=lut_uav, tau)

rmse<-function(x, y)  sqrt(mean((y-x)^2, na.rm=TRUE))

bias<-function(x, y)  mean(y-x, na.rm=TRUE)


# Load in the data --------------------------------------------------------

# Load the UAV and LiDAR images
uav_dsm<-raster("data/raster/uav_dsm_matched_cropped.tif")
uav_dtm<-raster("data/raster/uav_dtm_matched_cropped.tif")
lid_dsm<-raster("data/raster/lid_dsm_matched_cropped.tif")
lid_dtm<-raster("data/raster/lid_dtm_matched_cropped.tif")
lid_chm<-raster("data/raster/lid_chm_matched_cropped.tif")

lid_chm<-readAll(lid_chm)
lid_dtm<-readAll(lid_dtm)


# Ensure the UAV and LiDAR DEMs are vertically aligned (numbers from manual alignment):
uav_dsm<-uav_dsm-17.6
uav_dtm<-uav_dtm-17.6
# remove very low values from UAV DSM
uav_dsm[uav_dsm<0]<-NA
uav_dtm[uav_dtm<0]<-NA
# Calculates the canopy height models for the UAV
uav_chm<-uav_dsm-uav_dtm # UAV DTM


# Calculate biomass --------------------------------------------------------

lid_agb<-raster_agb(lid_chm, cover_h_thresh = 20, scale = 100)
uav_agb<-raster_agb(uav_chm, cover_h_thresh = 20, scale = 100)
uav_cor_agb<-raster_agb(uav_chm+4.14, cover_h_thresh = 20, scale = 100)


par(mfrow=c(2,2))
col.pal<-list(color = colorRampPalette(brewer.pal(9,"GnBu"))(10))$color
col.breaks<-seq(0, 80, length=length(col.pal)+1)
# LiDAR measured AGB
plot(lid_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
# UAV measured AGB
plot(uav_cor_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
col.pal<-list(color = colorRampPalette(brewer.pal(9,"RdBu"))(10))$color
col.breaks<-seq(-15, 15, length=length(col.pal)+1)
# Absoute difference
plot((lid_agb-uav_cor_agb), col = col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
col.pal<-list(color = colorRampPalette(brewer.pal(9,"RdBu"))(10))$color
col.breaks<-seq(-1, 1, length=length(col.pal)+1)
# Relative difference
plot((lid_agb-uav_cor_agb)/lid_agb, col = col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)


hist((lid_agb-uav_cor_agb))


plot(uav_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)

agb_grid <-
  data.frame(
    lid = values(lid_agb),
    uav = values(uav_cor_agb)
  )


ggplot(agb_grid, aes(x = lid, y = uav)) +
  geom_point() +
  #geom_smooth(method ="lm", se= FALSE, color = "darkgrey") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, color = "darkgrey") +
  geom_abline(intercept = 0, slope = 1)+
  xlim(0, 75) +
  ylim(0,75) +
  #xlab(expression(paste(AGB[LiDAR], " (Mg ", ha^-1, ")"))) +
  xlab(expression(paste("LiDAR AGB", " (Mg ", ha^-1, ")"))) +
  ylab(expression(paste("SFM AGB", " (Mg ", ha^-1, ")"))) +
  theme_minimal()



lid_max <-
  aggregate(
    lid_chm,
    fact = n_pix,
    fun = function(x, na.rm)
      max(x),
    progress = progress
  )

