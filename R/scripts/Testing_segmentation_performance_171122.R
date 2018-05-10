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
# Calculate the 90th percentile of the manual crowns height predicted from radius:
man_lid<-readOGR("data/shape/Manual trees LiDAR UTM.shp")
man_lid$Area<-as.numeric(as.character(man_lid$Area))
man_lid$R<-sqrt(man_lid$Area/pi) # Calculate a pseudoradius
man_lid$lid_CH_mean<-raster::extract(lid_chm, man_lid, fun=mean)
man_lid$lid_CH_max<-raster::extract(lid_chm, man_lid, fun=max)
lut_lid<-rq_lut(x=man_lid$lid_CH_max, y=man_lid$R, log=TRUE)

man_uav<-readOGR("data/shape/Manual trees UTM.shp")
man_uav$Area<-as.numeric(as.character(man_uav$Area))
man_uav$R<-sqrt(man_uav$Area/pi) # Calculate a pseudoradius
man_uav$uav_CH_mean<-raster::extract(uav_chm, man_uav, fun=mean)
man_uav$uav_CH_max<-raster::extract(uav_chm, man_uav, fun=max)
lut_uav<-rq_lut(x=man_uav$uav_CH_max, y=man_uav$R, log=TRUE)

# Comparison of manually segmented crown size:
png("Figures/AssessingForestRestoration/Radii.png",
    width=14,height=14, units="cm", res=300)
ggplot(data.frame(LiDAR = man_lid$R, SFM = man_uav$R), aes(x = LiDAR, y = SFM)) +
  geom_point() +
  geom_abline(linetype = "dashed") +
  xlab("LiDAR crown radius (m)") +
  ylab("SFM crown radius (m)") +
  theme_minimal()
dev.off()

rmse(x= man_lid$R, y = man_uav$R)
bias(x= man_lid$R, y = man_uav$R)
t.test(x= man_lid$R, y = man_uav$R, paired = TRUE)

# Look up tables to calculated height from radii:
lut_rtoh_lid<-rq_lut(x=man_lid$R, y=man_lid$lid_CH_max, log=TRUE)
lut_rtoh_uav<-rq_lut(x=man_uav$R, y=man_uav$uav_CH_max, log=TRUE)

# Load in the cost data for all of the parameters:
cost<-read.csv("data/shape/ITC_trees_params_cost_uav/Matched_costed/Cost.csv")


# Make a parameter matrix:
THRESHSeed_vec<-seq(from=0.3, to=0.9, by=0.1)
THRESHCrown_vec<-seq(from=0.3, to=0.9, by=0.1)
SOBELstr_vec<-seq(from=1.5, to=6, by=1)
params<-expand.grid(THRESHSeed=THRESHSeed_vec,
                    THRESHCrown=THRESHCrown_vec,
                    SOBELstr=SOBELstr_vec)


# Identifying best parameters ----------------------------------------------

plot(cost$mean.cost, type='l', ylim=c(0,1))
points(cost$med.cost, type='l', col='red')
points(cost$sd.cost, type='l', col='blue')

plot(cost$sd.cost~cost$med.cost, type='p', col='blue')

cost[(cost$med.cost < 0.2 & cost$sd.cost <0.25), ]
params[(cost$med.cost < 0.2 & cost$sd.cost <0.25), ]

par(mfrow=c(2,2), mar=c(4,4,2,2))
plot(cost[,'med.cost'], type='l'); abline(v = which.max(out2[,'n'] - out2[,'error']), col='blue')

cost[(cost$med.cost < 0.6 & cost$sd.cost <1), ]
params[(cost$med.cost < 0.6 & cost$sd.cost <1), ]





# Assessing performance of the best parameters ---------------------------

auto<-readOGR('data/shape/ITC_trees_params_cost_uav/seed_0.6_crown_0.7_sobel_5.5.shp')
man_crowns<-readOGR('data/shape/ITC_trees_params_cost_uav/Matched_costed/seed_0.6_crown_0.7_sobel_5.5.shp')

# Ensure cost is calcuated correctly
man_crowns$cost<-(man_crowns$overseg+man_crowns$undersg)/2 + man_crowns$tr_mtch
hist(man_crowns$cost); median(man_crowns$cost)

# Calculate shapebound of autocrowns
auto@data<-mutate(auto@data, shpbnd = hbnd_mx + hbnd_mn +soblbnd)

hist(auto$R) # Assess crown size distribution
conf_auto<-auto[auto$R>=0.5,] # Exclude tiny crowns
n_autocrowns<-nrow(conf_auto)

par(mfrow=c(1,2))
hist(crowns$allmbnd); hist(crowns$shpbnd)

ggplot(auto@data, aes(x=trHghts_mx, y=R, color=allmbnd)) +
  geom_point() +
  geom_point(data = allom_lid, color="red") +
  geom_point(data = allom_uav, color="blue") +
  geom_point(data = man_crowns@data, color=man_crowns@data$size_rt)

par(mfrow=c(1,2))
hist(man_crowns$overseg); median(man_crowns$overseg); sd(man_crowns$overseg)
hist(man_crowns$undersg); median(man_crowns$undersg); sd(man_crowns$undersg)

hist(man_crowns$size_rt); median(man_crowns$size_rt); sd(man_crowns$size_rt)

ggplot(man_crowns@data, aes(x=shpbnd, y=size_rt)) +geom_point()

shpbnd_thresh<-0.8
all_conf<-man_crowns@data %>% mutate(conf_group = "all")
high_conf<-man_crowns@data %>%
  filter(shpbnd >= shpbnd_thresh) %>%
  mutate(conf_group = "high")
man_crowns_by_conf<-rbind(all_conf, high_conf)
ggplot(man_crowns_by_conf, aes(size_rt, group=conf_group, color=conf_group)) +geom_density()
# This graph would be good to present.
man_crowns_by_conf %>%
  group_by(conf_group) %>%
  summarise(n=length(size_rt), mean(size_rt), mean(size_rt), median(size_rt), median(cost), sd(cost), sd(size_rt))

ggplot(man_crowns@data, aes(x=trHghts_mx, y=cost)) +
  geom_point() +
  geom_smooth(method="lm")
summary(lm(cost~trHghts_mx, data=man_crowns@data))

ggplot(man_crowns@data, aes(x=trHghts_mx, y=size_rt, size=shpbnd)) +
  geom_point() +
  geom_smooth(method="lm")


sum(man_crowns$shpbnd>=shpbnd_thresh)
hist(man_crowns$cost[man_crowns$shpbnd>=shpbnd_thresh]); median(man_crowns$cost[man_crowns$shpbnd>=shpbnd_thresh])
hist(man_crowns$size_rt[man_crowns$shpbnd>=shpbnd_thresh]); median(man_crowns$size_rt[man_crowns$shpbnd>=shpbnd_thresh])


man_crowns_centroid<-gCentroid(man_crowns, byid = TRUE, id = as.character(man_crowns$id))
man_crowns_coords<-coordinates(man_crowns_centroid)
xy_cost<-data.frame(x=man_crowns_coords[,'x'],
                    y=man_crowns_coords[,'y'],
                    z=man_crowns$cost)

xy_cost<-data.frame(x=man_crowns_coords[,'x'],
                    y=man_crowns_coords[,'y'],
                    z=man_crowns$size_rt)

fm.cost1 <-mgcv::gam(z ~ s(x, y, k=10),
                     data= xy_cost)
summary(fm.cost1)


# Assess autocrowns -------------------------------------------------------
conf_auto<-conf_auto[conf_auto$shpbnd>0.8,]
nrow(conf_auto)/n_autocrowns

# Plot out the segmented crowns:
par(mfrow=c(1,1), mar=c(0,0,0,0)); plot(uav_chm)
#uav_chm_crop<-raster::select(uav_chm)
e<-extent(313590.5, 313854.2, 9746456, 9746638)
uav_chm_crop<-crop(uav_chm, e)
col.pal<-list(color = colorRampPalette(brewer.pal(9,"GnBu"))(10))$color
col.breaks<-seq(0, 60, length=length(col.pal)+1)
png("Figures/AssessingForestRestoration/Figure1_itc.png",
    width=(e@xmax-e@xmin)/10, height=(e@ymax-e@ymin)/10, units="cm", res=500)
plot(uav_chm_crop, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, box=FALSE)
plot(auto, add=TRUE)
plot(conf_auto, add=TRUE, col = rgb(0.7, 0, 0.3, 0.5))
dev.off()

# Calculate spatial distribution
conf_auto_q<-mutate(conf_auto@data, qx = trunc(max_x/100), qy = trunc(max_y/100))
n_ha<-conf_auto_q %>% mutate(q = paste(qx, qy)) %>%
  group_by(q) %>%
  summarise(n=length(q)) %>%
  ungroup() %>%
  separate(q, into = c("x", "y"))

ggplot(n_ha, aes(x = x, y = y, size = n)) + geom_point() # spatial density
ggplot(n_ha, aes(n)) + geom_histogram() # frequency distribution
n_ha %>% summarise(mean(n, na.rm=TRUE), sd(n, na.rm=TRUE))





# Correct tree heights -------------------------------------

# Extract the lidar heights for the automatic trees
# This is likely to be slow: save the result and reload it:
conf_auto$trHghts_mx_lid<-raster::extract(lid_chm, conf_auto, fun=max)[,1]

# Is it appropriate to just exclude outliers? ... maybe not
conf_auto$trHghts_abs_diff<-conf_auto$trHghts_mx_lid - conf_auto$trHghts_mx
plot(lid_chm-uav_chm)
plot(conf_auto[conf_auto$trHghts_abs_diff>8 & conf_auto$trHghts_mx>20,], add = TRUE)

sum(is.na(conf_auto$trHghts_mx_lid))
# Exclude crowns without data for both
conf_auto<-conf_auto[!is.na(conf_auto$trHghts_mx_lid),]

# Allometric estimation
#conf_auto$trHghts_mx_uavallom50<-allom_lookup(conf_auto$R, lut_rtoh_uav, tau=50, antilog=TRUE)/2
#conf_auto$trHghts_mx_uavallom50<-rtoh_lookup(conf_auto$R, lut_uav, tau=60)
conf_auto$trHghts_mx_uavallom50<-rtoh_lookup(conf_auto$R, lut_lid, tau=50)


# Estimating the bias in the tree heights from the:
# ... CHMs
uav_chm10<-uav_chm
uav_chm10[uav_chm10<10]<-NA
lid_chm10<-lid_chm
lid_chm10[uav_chm<10]<-NA
bias(x = values(lid_chm10), y = values(uav_chm10))
# ... manually segmented crowns
bias_man<-bias(y=man_uav$uav_CH_max, x=man_lid$lid_CH_max)
# ... the automatically segmented crowns
bias_auto<-bias(y=conf_auto$trHghts_mx[conf_auto$trHghts_mx>10], x=conf_auto$trHghts_mx_lid[conf_auto$trHghts_mx>10])
bias_auto<-bias(y=conf_auto$trHghts_mx[conf_auto$trHghts_mx>10], x=conf_auto$trHghts_mx_lid[conf_auto$trHghts_mx>10])

plot(y=conf_auto$trHghts_mx[conf_auto$trHghts_mx>10], x=conf_auto$trHghts_mx_lid[conf_auto$trHghts_mx>10])
points(y=man_uav$uav_CH_max[man_uav$uav_CH_max>10], x=man_lid$lid_CH_max[man_uav$uav_CH_max>10], col = "red")
abline(0,1, col="blue")



# Estimating LiDAR height from SMF using a linear model
fm.hght2<-lm(trHghts_mx_lid~trHghts_mx, data=conf_auto)


fm.hght2_gam<-mgcv::gam(trHghts_mx_lid~s(trHghts_mx, k=2), data=conf_auto)
summary(fm.hght2_gam)
conf_auto$trHghts_mx_diff<-predict(fm.hght2_gam, newdata = data.frame(trHghts_mx=conf_auto$trHghts_mx))

conf_auto$trHghts_mx_diff<-conf_auto$trHghts_mx+5
conf_auto$trHghts_mx_bias_man<-conf_auto$trHghts_mx-bias_man
conf_auto$trHghts_mx_bias_auto<-conf_auto$trHghts_mx-bias_auto


# par(mfrow = c(1, 1))
# plot(
#   trHghts_mx_lid ~ trHghts_mx,
#   data = conf_auto,
#   xlim = c(0, 50),
#   ylim = c(0, 50)
# )
# abline(0, 1, col = 'blue')
# points(
#   x = conf_auto$trHghts_mx_lid,
#   y = predict(fm.hght2_gam),
#   col = 'red'
# )

# Adding the average paired difference between LiDAR and SFM height surfaces

# pc_diff <- conf_auto@data %>%
#   filter(trHghts_mx > 20) %>%
#   mutate(diff = trHghts_mx_lid - trHghts_mx) %>%
#   mutate(ratio = trHghts_mx / trHghts_mx_lid) %>%
#   summarise(diff = median(diff, na.rm = TRUE),
#             ratio = median(ratio, na.rm = TRUE))
#
# conf_auto$trHghts_mx_diff <- conf_auto$trHghts_mx * (1 / pc_diff$ratio)

g1<-ggplot(conf_auto@data, aes(x=trHghts_mx_lid, y=trHghts_mx)) +
  geom_point(alpha=0.04, size = 2) +
  #geom_smooth(method="lm", formula = y~poly(x, 2), color = "darkgrey", se = FALSE) +
  ylim(0, 45) +
  xlab("LiDAR height (m)") +
  #ylab(expression(paste(Height[SFM], " (m)"))) +
  ylab("SFM height (m)") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid") +
  theme_minimal()

fm_test1<-lm(trHghts_mx~trHghts_mx_lid, data = conf_auto@data)
fm_test2<-lm(trHghts_mx~poly(trHghts_mx_lid,2), data = conf_auto@data)
AIC(fm_test1, fm_test2)

g2<-ggplot(conf_auto@data, aes(x=trHghts_mx_lid, y=trHghts_mx_uavallom50)) +
  geom_point(alpha=0.03, size = 3) +
  geom_smooth(method="lm", color="darkgrey", se=FALSE) +
  ylim(0, 45) +
  xlab("") +
  ylab("") +
  #ylab(expression(paste(Height[SFMallom], " (m)"))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid") +
  theme_minimal()

g3<-ggplot(conf_auto@data, aes(x=trHghts_mx_lid, y=trHghts_mx_bias_man)) +
  geom_point(alpha=0.03, size = 3) +
  geom_smooth(method="lm", color="darkgrey", se=FALSE) +
  ylim(0, 45) +
  #xlab(expression(paste(Height[LiDAR], " (m)"))) +
  xlab("") +
  ylab("") +
  #ylab(expression(paste(Height[SFM+alpha], " (m)"))) +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "solid") +
  theme_minimal()

tree_height_rmse_lid_uav<-boot_height_rmse(x=conf_auto$trHghts_mx_lid, y=conf_auto$trHghts_mx, height_window = 5,
                                           n_boot = 1000,
                                           n_sample = 30)

tree_height_rmse_lid_uavallom50<-boot_height_rmse(x=conf_auto$trHghts_mx_lid, y=conf_auto$trHghts_mx_uavallom50, height_window = 5,
                                           n_boot = 1000,
                                           n_sample = 30)

tree_height_rmse_lid_uav_bias_man<-boot_height_rmse(x=conf_auto$trHghts_mx_lid, y=conf_auto$trHghts_mx_bias_man, height_window = 5,
                                           n_boot = 1000,
                                           n_sample = 30)



g1_rmse<-ggplot(tree_height_rmse_lid_uav, aes(x = height, y = (rmse/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_rmse/100)*height, ymax = (upper_rmse/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  #xlim(0, 32) +
  #ylim(0, 75) +
  xlab("LiDAR height (m)") +
  #ylab(expression(paste("% ", RMSE[SFM]))) +
  ylab("RMSE (m)") +
  theme_minimal()


g1_bias<-ggplot(tree_height_rmse_lid_uav, aes(x = height, y = (bias/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_bias/100)*height, ymax = (upper_bias/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  #xlim(0, 32) +
  #ylim(0, 75) +
  xlab("LiDAR height (m)") +
  ylab("Bias (m)") +
  theme_minimal()

g2_rmse<-ggplot(tree_height_rmse_lid_uavallom50, aes(x = height, y = (rmse/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_rmse/100)*height, ymax = (upper_rmse/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "dashed", color = "red") +
  #xlim(0, 32) +
  #ylim(0, 75) +
  xlab("LiDAR height (m)") +
  #ylab(expression(paste("% ", RMSE[SFM]))) +
  ylab("") +
  theme_minimal()


g2_bias<-ggplot(tree_height_rmse_lid_uavallom50, aes(x = height, y = (bias/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_bias/100)*height, ymax = (upper_bias/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "dashed", color = "red") +
  #xlim(0, 32) +
  #ylim(0, 75) +
  xlab("Height (m)") +
  ylab("") +
  theme_minimal()

g3_rmse<-ggplot(tree_height_rmse_lid_uav_bias_man, aes(x = height, y = (rmse/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_rmse/100)*height, ymax = (upper_rmse/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "dashed") +
  #xlim(0, 32) +
  #ylim(0, 75) +
  xlab("") +
  #ylab(expression(paste("% ", RMSE[SFM]))) +
  ylab("") +
  theme_minimal()


g3_bias<-ggplot(tree_height_rmse_lid_uav_bias_man, aes(x = height, y = (bias/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_bias/100)*height, ymax = (upper_bias/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "dashed") +
  #xlim(0, 32) +
  #ylim(0, 75) +
  xlab("Height (m)") +
  ylab("") +
  theme_minimal()

#lay<-matrix(1:9, ncol=3, byrow = TRUE)
png("Figures/AssessingForestRestoration/Figure2_itcHeights_b.png",
    width=25,height=7.5, units="cm", res=500)
grid.arrange(g1, g1_rmse, g1_bias, ncol=3)
dev.off()

# png("Figures/AssessingForestRestoration/Figure2_itcHeights.png",
#     width=20,height=20, units="cm", res=500)
# grid.arrange(g1, g2, g3, g1_rmse, g2_rmse, g3_rmse, g1_bias, g2_bias, g3_bias, ncol=3)
# dev.off()


cor.test(y=conf_auto@data$trHghts_mx, x=conf_auto@data$trHghts_mx_lid)
cor.test(y=conf_auto@data$trHghts_mx_uavallom50, x=conf_auto@data$trHghts_mx_lid)
cor.test(y=conf_auto@data$trHghts_mx_diff, x=conf_auto@data$trHghts_mx_lid)

rmse(y=conf_auto@data$trHghts_mx, x=conf_auto@data$trHghts_mx_lid)
rmse(y=conf_auto@data$trHghts_mx_uavallom50, x=conf_auto@data$trHghts_mx_lid)
rmse(y=conf_auto@data$trHghts_mx_bias_man, x=conf_auto@data$trHghts_mx_lid)

bias(y=conf_auto@data$trHghts_mx, x=conf_auto@data$trHghts_mx_lid)
bias(y=conf_auto@data$trHghts_mx_uavallom50, x=conf_auto@data$trHghts_mx_lid)
bias(y=conf_auto@data$trHghts_mx_bias_man, x=conf_auto@data$trHghts_mx_lid)

allom_uav90<-data.frame(trHghts_mx=0:50, R=(htod_lookup(0:50, lut=lut_uav, tau=90)/2))
allom_lid90<-data.frame(trHghts_mx=0:50, R=(htod_lookup(0:50, lut=lut_lid, tau=90)/2))
allom_lid50<-data.frame(trHghts_mx=0:50, R=(htod_lookup(0:50, lut=lut_lid, tau=50)/2))
allom_uavR_allomH<-dplyr::select(conf_auto@data, trHghts_mx=trHghts_mx_uavallom50, R=R)
allom_uavR_lidH<-dplyr::select(conf_auto@data, trHghts_mx=trHghts_mx_lid, R=R)
allom_uavR_diff<-dplyr::select(conf_auto@data, trHghts_mx=trHghts_mx_diff, R=R)

ggplot(auto@data, aes(x=trHghts_mx, y=R)) +
  geom_point(alpha = 0.5) +
  geom_point(data = allom_uavR_allomH, colour = "orange", alpha = 0.5) +
  geom_point(data = allom_uavR_diff, colour = "yellow", alpha = 0.5) +
  #geom_point(data = allom_uavR_lidH, colour = "green", alpha = 0.5) +
  geom_point(data = allom_uav90, color="blue") +
  geom_point(data = allom_lid90, color="red") +
  geom_point(data = allom_lid50, color="red")


# Predict DTM from the estimated tree heights ----

k=1000 # A hack as the k object doesn't seem to be able to be read as a functional argument

# Allometric correction
est_dtm_allom<-dtm_predict(x=conf_auto$max_x, y=conf_auto$max_y, z=conf_auto$trHghts_mx_uavallom50,
                     img=uav_dsm,
                     predict_grid_by=25,
                     k=k)

# Difference correction
est_dtm_diff<-dtm_predict(x=conf_auto$max_x, y=conf_auto$max_y, z=conf_auto$trHghts_mx_bias_man,
                     img=uav_dsm,
                     predict_grid_by=25,
                     k=k)

# Only lowering the DTM / not allowing increases:
#low.est_dtm_allom<-uav_dtm
#low.est_dtm_allom[low.est_dtm_allom>est_dtm]<-est_dtm_allom[low.est_dtm_allom>est_dtm_allom] # Calculate the low estimate CHM

# Plot the DTMs ----

par(mfrow=c(1,4))
col.pal<-rev(list(color = colorRampPalette(brewer.pal(9,"YlOrBr"))(10))$color)
col.breaks<-seq(-1, 60, length=length(col.pal)+1)
plot(lid_dtm, col=col.pal, breaks=col.breaks, colNA='black')
plot(uav_dtm, col=col.pal, breaks=col.breaks, colNA='black')
plot(est_dtm_allom, col=col.pal, breaks=col.breaks, colNA='black')
plot(est_dtm_diff, col=col.pal, breaks=col.breaks, colNA='black')

# Caculate the CHMs ----
uav_est_allom_chm<-uav_dsm-est_dtm_allom # Calculate the estimated CHM
uav_est_diff_chm<-uav_dsm-est_dtm_diff # Calculate the estimated CHM
#uav_low.est_allom_chm<-uav_dsm-low.est_dtm # Calculate the estimated CHM

# Do not permit values to drop beneath 0
uav_chm[uav_chm<0]<-0
uav_est_allom_chm[uav_est_allom_chm<0]<-0
uav_est_diff_chm[uav_est_diff_chm<0]<-0
#uav_low.est_allom_chm[uav_low.est_allom_chm<0]<-0


mean(values(lid_chm), na.rm=TRUE); sd(values(lid_chm), na.rm=TRUE)
mean(values(uav_chm), na.rm=TRUE); sd(values(uav_chm), na.rm=TRUE)
mean(values(uav_est_allom_chm), na.rm=TRUE); sd(values(uav_est_allom_chm), na.rm=TRUE)
mean(values(uav_est_diff_chm), na.rm=TRUE); sd(values(uav_est_diff_chm), na.rm=TRUE)

grid_scale <- 10

lid_chm_grid<-grid_fun(lid_chm, grid_by=grid_scale, fun=function(x) mean(x, na.rm=TRUE))
uav_chm_grid<-grid_fun(uav_chm, grid_by=grid_scale, fun=function(x) mean(x, na.rm=TRUE))
uav_est_allom_chm_grid<-grid_fun(uav_est_allom_chm, grid_by=grid_scale, fun=function(x) mean(x, na.rm=TRUE))
uav_est_diff_chm_grid<-grid_fun(uav_est_diff_chm, grid_by=grid_scale, fun=function(x) mean(x, na.rm=TRUE))

chm_grid <-
  data.frame(
    lid = values(lid_chm_grid),
    uav = values(uav_chm_grid),
    uav_est_allom = values(uav_est_allom_chm_grid),
    uav_est_diff = values(uav_est_diff_chm_grid)
  )


# Produce height density figure ----

# !!!! BE WARNED - THIS TAKES A FEW MINUTES !!!!

m1<-data.frame(model = "lid_chm", height = values(lid_chm))
m2<-data.frame(model = "uav_chm", height = values(uav_chm))
m3<-data.frame(model = "uav_est_allom_chm", height = values(uav_est_allom_chm))
m4<-data.frame(model = "uav_est_diff_chm", height = values(uav_est_diff_chm))
m1234<-rbind(m1, m2, m3, m4)
m1234<-m1234[complete.cases(m1234),]
rm(m1);rm(m2);rm(m3);rm(m4);

library(palettetown)

mypal<-"arcanine" %>% ichooseyou(4)
mypal<-"bulbasaur" %>% ichooseyou(4)


mypal<-list(color = brewer.pal(4, "Set1"))$color

png("Figures/AssessingForestRestoration/Figure5.png",
    width=15,height=10, units="cm", res=500)
ggplot(m1234, aes(height, group=model, color=model)) +
  geom_density(alpha = 0.3) +
  #stat_density(geom="line") +
  xlab("Height") +
  ylab("Density") +
  scale_color_manual( values = mypal) +
  theme_minimal()
dev.off()

rmse(chm_grid$uav, chm_grid$lid)
rmse(chm_grid$uav_est_allom, chm_grid$lid)
rmse(chm_grid$uav_est_diff, chm_grid$lid)

bias(chm_grid$uav, chm_grid$lid)
bias(chm_grid$uav_est_allom, chm_grid$lid)
bias(chm_grid$uav_est_diff, chm_grid$lid)



g5<-ggplot(chm_grid, aes(x=lid, y=uav)) +
  geom_point(alpha=0.02, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red") +
  geom_smooth(method = "lm", formula = y~poly(x, 3), color = "darkgrey", se = FALSE) +
  xlim(0, 32) +
  ylim(0, 35) +
  #ylab(expression(paste(Height[SFM], " (m)"))) +
  xlab("") +
  ylab("SFM Height (m)") +
  theme_minimal()

g6<-ggplot(chm_grid, aes(x=lid, y=uav_est_allom)) +
  geom_point(alpha=0.02, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red") +
  geom_smooth(method = "lm", formula = y~poly(x, 3), color = "darkgrey", se = FALSE) +  xlim(0, 32) +
  ylim(0, 35) +
  #ylab(expression(paste(Height[SFMallom], " (m)"))) +
  xlab("") +
  ylab("") +
  theme_minimal()

g7<-ggplot(chm_grid, aes(x=lid, y=uav_est_diff)) +
  geom_point(alpha=0.02, size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "red") +
  geom_smooth(method = "lm", formula = y~poly(x, 3), color = "darkgrey", se = FALSE) +  xlim(0, 32) +
  ylim(0, 35) +
  #ylab(expression(paste(Height[SFM+alpha], " (m)"))) +
  xlab("") +
  ylab("") +
  theme_minimal()

# Bootstrap the RMSE per height ----

height_rmse_lid_uav<-boot_height_rmse(x=chm_grid$lid, y=chm_grid$uav, height_window = 5,
            n_boot = 1000,
            n_sample = 30)

height_rmse_lid_uav_est_allom<-boot_height_rmse(x=chm_grid$lid, y=chm_grid$uav_est_allom, height_window = 5,
                                      n_boot = 1000,
                                      n_sample = 30)

height_rmse_lid_uav_est_diff<-boot_height_rmse(x=chm_grid$lid, y=chm_grid$uav_est_diff, height_window = 5,
                                                n_boot = 1000,
                                                n_sample = 30)


g_rmse1<-ggplot(height_rmse_lid_uav, aes(x = height, y = (rmse/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_rmse/100)*height, ymax = (upper_rmse/100)*height), alpha = 0.4) +
  #geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  xlim(0, 32) +
  ylim(0, 12) +
  xlab("") +
  #ylab(expression(paste("% ", RMSE[SFM]))) +
  ylab("RMSE (m)") +
  theme_minimal()
g_rmse2<-ggplot(height_rmse_lid_uav_est_allom, aes(x = height, y = (rmse/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_rmse/100)*height, ymax = (upper_rmse/100)*height), alpha = 0.4) +
  #geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  xlim(0, 32) +
  ylim(0, 12) +
  xlab("") +
  ylab("") +
  #ylab(expression(paste("% ", RMSE[SFMallom]))) +
  theme_minimal()
g_rmse3<-ggplot(height_rmse_lid_uav_est_diff, aes(x = height, y = (rmse/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_rmse/100)*height, ymax = (upper_rmse/100)*height), alpha = 0.4) +
  #geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  xlim(0, 32) +
  ylim(0, 12) +
  xlab("") +
  #ylab(expression(paste("% ", RMSE[SFM+alpha]))) +
  ylab("") +
  theme_minimal()


g_bias1<-ggplot(height_rmse_lid_uav, aes(x = height, y = (bias/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_bias/100)*height, ymax = (upper_bias/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  xlim(0, 32) +
  ylim(-10, 5) +
  #xlab(expression(paste(Height[LiDAR], " (m)"))) +
  xlab("LiDAR height (m)") +
  #ylab(expression(paste("% ", Bias[SFM]))) +
  ylab("Bias (m)") +
  theme_minimal()
g_bias2<-ggplot(height_rmse_lid_uav_est_allom, aes(x = height, y = (bias/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_bias/100)*height, ymax = (upper_bias/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  xlim(0, 32) +
  ylim(-10, 5) +
  #xlab(expression(paste(Height[LiDAR], " (m)"))) +
  xlab("LiDAR height (m)") +
  #ylab(expression(paste("% ", Bias[SFMallom]))) +
  ylab("") +
  theme_minimal()
g_bias3<-ggplot(height_rmse_lid_uav_est_diff, aes(x = height, y = (bias/100)*height)) +
  geom_line() +
  geom_ribbon(aes(ymin = (lower_bias/100)*height, ymax = (upper_bias/100)*height), alpha = 0.4) +
  geom_abline(intercept = 0, slope = 0, linetype= "solid", color = "red") +
  xlim(0, 32) +
  ylim(-10,5) +
  #xlab(expression(paste(Height[LiDAR], " (m)"))) +
  xlab("LiDAR height (m)") +
  #ylab(expression(paste("% ", Bias[SFM+alpha]))) +
  ylab("") +
  theme_minimal()

# Plot summary figures of CHM surfaces with RMSE & bias ----

png("Figures/AssessingForestRestoration/Figure3.png",
    width=20,height=20, units="cm", res=500)

grid.arrange(g5, g7, g6,
             g_rmse1, g_rmse3, g_rmse2,
             g_bias1, g_bias3, g_bias2, ncol=3)
dev.off()


# CHMs and their differences ----

chm_diff<-uav_chm_grid-lid_chm_grid
chm_est_allom<-uav_est_allom_chm_grid-lid_chm_grid
chm_est_diff<-uav_est_diff_chm_grid-lid_chm_grid

chm_grid$diff<-values(chm_diff)
chm_grid$allom_diff<-values(chm_est_allom)
chm_grid$diff_diff<-values(chm_est_diff)

chm_grid %>% dplyr::select(diff, allom_diff, diff_diff) %>%
  summarise(diff = quantile(diff, probs =0.95, na.rm = TRUE),
                allom_diff = quantile(allom_diff, probs =0.95, na.rm = TRUE),
                diff_diff = quantile(diff_diff, probs =0.95, na.rm = TRUE))

png("Figures/AssessingForestRestoration/Figure4.png",
    width=28,height=20, units="cm", res=500)
par(mfrow=c(2,4), mar=c(0,0,0,0), oma=c(0,0,0,0))
col.pal<-list(color = colorRampPalette(brewer.pal(9,"GnBu"))(10))$color
col.breaks<-seq(0, 60, length=length(col.pal)+1)
plot(lid_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
plot(uav_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
plot(uav_est_allom_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
plot(uav_est_diff_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
col.pal<-list(color = colorRampPalette(brewer.pal(11,"RdYlBu"))(10))$color
col.breaks<-seq(-10, 10, length=length(col.pal)+1)
plot(0, type="n", axes=FALSE, ylab="", xlab="")
plot(chm_diff, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
plot(chm_est_allom, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
plot(chm_est_diff, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=FALSE, box=FALSE)
dev.off()


# To produce the legend
par(mfrow=c(1,2), mar=c(0,0,0,0), oma=c(0,0,0,0))
col.pal<-list(color = colorRampPalette(brewer.pal(9,"GnBu"))(10))$color
col.breaks<-seq(0, 60, length=length(col.pal)+1)
plot(lid_chm, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
col.pal<-list(color = colorRampPalette(brewer.pal(11,"RdYlBu"))(10))$color
col.breaks<-seq(-10, 10, length=length(col.pal)+1)
plot(chm_diff, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)



# Spatial autocorrelation ----
Moran(chm_diff)
Moran(chm_est_allom)
Moran(chm_est_diff)

# Biomass ----

# Functions for calculating AGB taken from Jucker et al., 2017
# "A regional model for estimating the aboveground carbon density
# of Borneo's tropical forests from airborne laser scanning"

lid_agb<-raster_agb(lid_chm, cover_h_thresh = 20, scale = 100)
uav_agb<-raster_agb(uav_chm, cover_h_thresh = 20, scale = 100)
uav_est_allom_agb<-raster_agb(uav_est_allom_chm, cover_h_thresh = 20, scale = 100)
uav_est_diff_agb<-raster_agb(uav_est_diff_chm, cover_h_thresh = 20, scale = 100)


par(mfrow=c(1,4))
col.pal<-list(color = colorRampPalette(brewer.pal(9,"GnBu"))(10))$color
col.breaks<-seq(0, 80, length=length(col.pal)+1)
plot(lid_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
plot(uav_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
plot(uav_est_allom_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
plot(uav_est_diff_agb, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)

agb_grid <-
  data.frame(
    lid = values(lid_agb),
    uav = values(uav_agb),
    uav_est_allom = values(uav_est_allom_agb),
    uav_est_diff = values(uav_est_diff_agb)
  )

g_agb1<-ggplot(agb_grid, aes(x = lid, y = uav)) +
  geom_point() +
  #geom_smooth(method ="lm", se= FALSE, color = "darkgrey") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, color = "darkgrey") +
  geom_abline(intercept = 0, slope = 1)+
  xlim(0, 80) +
  ylim(0,125) +
  #xlab(expression(paste(AGB[LiDAR], " (Mg ", ha^-1, ")"))) +
  xlab(expression(paste("LiDAR AGB", " (Mg ", ha^-1, ")"))) +
  ylab(expression(paste("SFM AGB", " (Mg ", ha^-1, ")"))) +
  theme_minimal()

g_agb2<-ggplot(agb_grid, aes(x = lid, y = uav_est_allom)) +
  geom_point() +
  #geom_smooth(method ="lm", se= FALSE, color = "darkgrey") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, color = "darkgrey") +
  geom_abline(intercept = 0, slope = 1)+
  xlim(0, 80) +
  ylim(0,125) +
  xlab(expression(paste("LiDAR AGB", " (Mg ", ha^-1, ")"))) +
  ylab(expression(paste("SFM AGB", " (Mg ", ha^-1, ")"))) +
  theme_minimal()

g_agb3<-ggplot(agb_grid, aes(x = lid, y = uav_est_diff)) +
  geom_point() +
  #geom_smooth(method ="lm", se= FALSE, color = "darkgrey") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, color = "darkgrey") +
  geom_abline(intercept = 0, slope = 1)+
  xlim(0, 80) +
  ylim(0,125) +
  xlab(expression(paste("LiDAR AGB", " (Mg ", ha^-1, ")"))) +
  ylab(expression(paste("SFM AGB", " (Mg ", ha^-1, ")"))) +
  theme_minimal()


png("Figures/AssessingForestRestoration/Figure6.png",
    width=20,height=10, units="cm", res=500)
grid.arrange(g_agb1, g_agb3, g_agb2, ncol = 3)
dev.off()

rmse(agb_grid$lid, agb_grid$uav) / mean(agb_grid$lid, na.rm = TRUE)
rmse(agb_grid$lid, agb_grid$uav_est_allom) / mean(agb_grid$lid, na.rm = TRUE)
rmse(agb_grid$lid, agb_grid$uav_est_diff) / mean(agb_grid$lid, na.rm = TRUE)


bias(agb_grid$lid, agb_grid$uav) / mean(agb_grid$lid, na.rm = TRUE)
bias(agb_grid$lid, agb_grid$uav_est_allom) / mean(agb_grid$lid, na.rm = TRUE)
(bias(agb_grid$lid, agb_grid$uav_est_diff) / mean(agb_grid$lid, na.rm = TRUE)) *100

# Topographical position ----

# Calculate TPI from the LiDAR DTM:
#lid_tpi<-tpi(lid_dtm, scale = 75)
lid_dtm_grid<-grid_fun(lid_dtm, grid_by=grid_scale, fun=function(x) mean(x, na.rm=TRUE))
#lid_tpi_grid<-grid_fun(lid_tpi, grid_by=grid_scale, fun=function(x) mean(x, na.rm=TRUE))

lid_tpi_grid<-tpi(lid_dtm_grid, scale = 15)

plot(lid_dtm_grid)
plot(lid_tpi_grid)

lid_chm_grid_no.small<-lid_chm_grid
lid_chm_grid_no.small[lid_chm_grid_no.small<2]<-NA
chm_reldiff<-chm_diff/lid_chm_grid_no.small


par(mfrow=c(1,3), mar=c(4,4,2,2))
col.pal<-list(color = colorRampPalette(brewer.pal(11,"RdYlBu"))(10))$color
col.breaks<-seq(-10, 10, length=length(col.pal)+1)
plot(lid_tpi_grid, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
col.pal<-list(color = colorRampPalette(brewer.pal(11,"RdYlBu"))(10))$color
col.breaks<-seq(-30, 15, length=length(col.pal)+1)
plot(chm_diff, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)
col.pal<-list(color = colorRampPalette(brewer.pal(11,"RdYlBu"))(10))$color
col.breaks<-seq(-1, 1, length=length(col.pal)+1)
plot(chm_reldiff, col=col.pal, breaks=col.breaks, colNA='black', axes=FALSE, legend=TRUE, box=FALSE)


chm_grid$tpi = values(lid_tpi_grid)
chm_grid<-cbind(chm_grid, expand.grid(x=1:ncol(lid_chm_grid), y = 1:nrow(lid_chm_grid)))
chm_grid$uav.std<-scale(chm_grid$uav)[,1]
chm_grid$lid.std<-scale(chm_grid$lid)[,1]
chm_grid$tpi.std<-scale(chm_grid$tpi)[,1]

chm_grid_complete<-chm_grid[complete.cases(chm_grid),]

fm.mgcv1<-mgcv::gam(diff_diff~poly(tpi,2) + s(lid,k=3) + s(x, y, k=200), data = chm_grid_complete)
#fm.mgcv2<-mgcv::gam(uav.std~lid.std + tpi.std + s(x, y, k=1000), data = chm_grid)
#fm.mgcv3<-mgcv::gam(uav.std~lid.std + s(x, y, k=1000), data = chm_grid)
summary(fm.mgcv1)
anova(fm.mgcv1)
AIC(fm.mgcv1, fm.mgcv2, fm.mgcv3)


g_tpi1<-chm_grid %>%
  filter(lid>=2) %>%
  ggplot(aes(x=tpi, y=diff)) +
  geom_point(alpha = 0.03, size = 3) +
#  geom_smooth(method="lm") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, col = "darkgrey") +
  geom_abline(intercept = 0, slope = 0, col = "darkgrey", linetype = "dashed") +
  ylim(-20, 20) +
  xlab("Topographic position index (AU)") +
  ylab("Error (m)") +
  theme_minimal()

g_tpi2<-chm_grid %>%
  filter(lid>=2) %>%
  ggplot(aes(x=tpi, y=allom_diff)) +
  geom_point(alpha = 0.03, size = 3) +
  #  geom_smooth(method="lm") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, col = "darkgrey") +
  geom_abline(intercept = 0, slope = 0, col = "darkgrey", linetype = "dashed") +
  ylim(-20, 20) +
  xlab("Topographic position index (AU)") +
  ylab("") +
  theme_minimal()

g_tpi3<-chm_grid %>%
  filter(lid>=2) %>%
  ggplot(aes(x=tpi, y=diff_diff)) +
  geom_point(alpha = 0.03, size = 3) +
  #  geom_smooth(method="lm") +
  stat_smooth(method = "lm", se = FALSE, formula = y ~ poly(x, 2), size = 1, col = "darkgrey") +
  geom_abline(intercept = 0, slope = 0, col = "darkgrey", linetype = "dashed") +
  ylim(-20, 20) +
  xlab("Topographic position index (AU)") +
  ylab("") +
  theme_minimal()

png("Figures/AssessingForestRestoration/Figure7.png",
    width=20,height=8, units="cm", res=500)
grid.arrange(g_tpi1, g_tpi3, g_tpi2, ncol = 3)
dev.off()

fm.tpi1<-lm(diff~tpi, data = chm_grid_complete)
fm.tpi2<-lm(diff~poly(tpi, 2), data = chm_grid_complete)
AIC(fm.tpi1, fm.tpi2)
summary(fm.tpi1)
summary(fm.tpi2)

fm.tpi_allom1<-lm(allom_diff~tpi, data = chm_grid_complete)
fm.tpi_allom2<-lm(allom_diff~poly(tpi, 2), data = chm_grid_complete)
AIC(fm.tpi_allom1, fm.tpi_allom2)
summary(fm.tpi_allom1)
summary(fm.tpi_allom2)

fm.tpi_diff1<-lm(diff_diff~tpi, data = chm_grid_complete)
fm.tpi_diff2<-lm(diff_diff~poly(tpi, 2), data = chm_grid_complete)
AIC(fm.tpi_diff1, fm.tpi_diff2)
summary(fm.tpi_diff1)
summary(fm.tpi_diff2)




