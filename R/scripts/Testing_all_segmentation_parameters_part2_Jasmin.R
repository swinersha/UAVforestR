# Assing overlap with each tree -----------------
#
# Finding the trees with the greatest overlap compared
# with a set of manually segemented crowns
#
# For each set of parameters a shapefile  of the crowns that match
# is saved of the segmented trees

rm(list=ls())

library(raster)
library(maptools)
library(rgdal)
library(ggplot2)
library(rgeos)
library(foreach)
library(doParallel)

# Load R source files
R_source_files<-list.files(
  path = "R",
  pattern = "*.R$",
  full.names = TRUE
)
sapply(R_source_files, function(x) source(x, local = FALSE,  echo = FALSE))

lid_chm<-raster("data/raster/lid_chm_matched_cropped.tif")
manual<-readOGR("data/shape/Manual trees UTM.shp")
crs(manual)<-crs(lid_chm)

# out<-matrix(0, nrow=nrow(params), ncol=9)
# colnames(out)<-c('manual.id',
#                  'n.unmatch',
#                  'mean.over',
#                  'sd.over',
#                  'mean.under',
#                  'sd.under',
#                  'mean.cost',
#                  'med.cost',
#                  'sd.cost')

# Make a parameter matrix:
THRESHSeed_vec<-seq(from=0.3, to=0.9, by=0.1)
THRESHCrown_vec<-seq(from=0.3, to=0.9, by=0.1)
SOBELstr_vec<-seq(from=1.5, to=6, by=1)
params<-expand.grid(THRESHSeed=THRESHSeed_vec,
                    THRESHCrown=THRESHCrown_vec,
                    SOBELstr=SOBELstr_vec)

cl<-parallel::makeCluster(2)
registerDoParallel(cl)

# cost<-foreach(i = 1:nrow(params),
#         .packages = c("raster", "rgeos", "sp", "rgdal"),
#         .combine = rbind,
#         .inorder = TRUE) %dopar% {

out<-vector(mode="list", length=nrow(params))
for(i in 1:nrow(params)){
  file_name<-paste('data/shape/ITC_trees_params_cost_uav/seed_', params[i,1],
                   '_crown_', params[i,2],
                   '_sobel_', params[i,3], '.shp', sep='')
  print(file_name)
  ut<-readOGR(file_name)
  crs(ut)<-crs(lid_chm)
  # Find the crowns with the greatest overlap:
  matched<-crown_overlap(auto_trees=ut, manual_trees=manual, buffer_by=60)

  # Summarises the quality of the match:
  # out[i,1]<-i
  # out[i,2]<-sum(matched$tree_match)
  # out[i,3]<-mean(matched$overseg)
  # out[i,4]<-sd(matched$overseg)
  # out[i,5]<-mean(matched$underseg)
  # out[i,6]<-sd(matched$underseg)
  # out[i,7]<-mean(matched$cost)
  # out[i,8]<-median(matched$cost)
  # out[i,9]<-sd(matched$cost)

  out_tmp<-data.frame(manual.id=i,
                  n.unmatch = sum(matched$tree_match),
                  mean.over = mean(matched$overseg),
                  sd.over = sd(matched$overseg),
                  mean.under = mean(matched$underseg),
                  sd.under = sd(matched$underseg),
                  mean.cost = mean(matched$cost),
                  med.cost = median(matched$cost),
                  sd.cost = sd(matched$cost),
                  mean.size_ratio = mean(matched$size_ratio),
                  med.size_ratio = median(matched$size_ratio),
                  sd.size_ratio = sd(matched$size_ratio)
                  )

  out[[i]]<-out_tmp
  # add the auto tree data to the manual trees:
  matched<-matched[!is.na(matched$id_auto_tree),]
  auto<-ut[matched$id_auto_tree,]
  auto@data$shpbnd<-rowSums(auto@data[,c('hbnd_mn', 'hbnd_mx', 'soblbnd')])
  auto@data$spcbnd<-rowSums(auto@data[,c('allmbnd', 'crwnbnd')])
  auto@data<-auto@data[,c('trHghts_mn', 'trHghts_mx', 'R', 'shpbnd', 'spcbnd', 'sobl_mn')]
  matched@data<-cbind(matched@data, auto@data)

  # par(mfrow=c(1,1), mar=c(0,0,0,0))
  # plot(matched)
  # plot(auto, add=TRUE, border='red')

  # Save the output:
  file_name<-paste('data/shape/ITC_trees_params_cost_uav/Matched_costed/seed_', params[i,1],
                   '_crown_', params[i,2],
                   '_sobel_', params[i,3], sep='')
  writeOGR(matched,
           dsn = paste(file_name, '.shp', sep=''),
           layer = basename(file_name),
           drive = 'ESRI Shapefile')

  #return(out)
}
#stopCluster(cl)
cost<-do.call(rbind, out)
cost<-as.data.frame(cost)

write.csv(cost, "data/shape/ITC_trees_params_cost_uav/Matched_costed/Cost.csv")

cost<-read.csv("data/shape/ITC_trees_params_cost_uav/Matched_costed/Cost.csv")

params_cost<-cbind(params, cost)
params_cost$cov.cost<-params_cost$sd.cost/params_cost$mean.cost


pdf("Figures/AssessingForestRestoration/Appendix_Figure1_sd.pdf",
    width=6,height=6)
ggplot(params_cost, aes(x = THRESHSeed, y = med.cost, color = as.factor(SOBELstr))) +
  geom_point() +
  geom_smooth(method = "loess", se = FALSE, size = 0.3) +
  #ylim(0.15, 0.45) +
  scale_color_grey(start = 0.8, end = 0.2) +
  facet_wrap(~as.factor(THRESHCrown), ncol = 2) +
  theme_minimal()
dev.off()

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

points(out2[,'mean.overlap'], type='l', col='red')
points(out2[,'error'], type='l', col='green')



par(mfrow=c(1,2), mar=c(4,4,2,2))
plot(out[,'mean.ratio'], type='l')
points(out[,'mean.overlap'], type='l', col='red')
plot(out[,'sd.ratio'], type='l')
points(out[,'sd.overlap'], type='l', col='red')

crowns<-readOGR('data/shape/ITC_trees_params_cost_uav/seed_0.6_crown_0.7_sobel_5.5.shp')
man_crowns<-readOGR('data/shape/ITC_trees_params_cost_uav/Matched_costed/seed_0.6_crown_0.7_sobel_5.5.shp')

man_crowns$cost<-(man_crowns$overseg+man_crowns$undersg)/2 + man_crowns$tr_mtch

hist(man_crowns$cost); median(man_crowns$cost)

allom_lid<-data.frame(trHghts_mx=0:50, R=(htod(0:50, tau=90)/2))
allom_uav<-data.frame(trHghts_mx=0:50, R=(htod_lookup(0:50, lut=lut_uav, tau=90)/2))

ggplot(crowns@data, aes(x=trHghts_mx, y=R)) +
  geom_point() +
  geom_point(data = allom_lid, color="red") +
  geom_point(data = allom_uav, color="blue")

hist(crowns$allmbnd)

ggplot(crowns@data, aes(x=trHghts_mx, y=R, color=allmbnd)) +
  geom_point() +
  geom_point(data = allom_lid, color="red") +
  geom_point(data = allom_uav, color="blue") +
  geom_point(data = man_crowns@data, color=man_crowns@data$size_rt)

par(mfrow=c(1,2))
hist(man_crowns$overseg); median(man_crowns$overseg); sd(man_crowns$overseg)
hist(man_crowns$undersg); median(man_crowns$undersg); sd(man_crowns$undersg)

hist(man_crowns$size_rt); median(man_crowns$size_rt); sd(man_crowns$size_rt)

ggplot(man_crowns@data, aes(x=shpbnd, y=size_rt)) +geom_point()

shpbnd_thresh<-1
all_conf<-man_crowns@data %>% mutate(conf_group = "all")
high_conf<-man_crowns@data %>%
  filter(shpbnd >= shpbnd_thresh) %>%
  mutate(conf_group = "high")
man_crowns_by_conf<-rbind(all_conf, high_conf)
ggplot(man_crowns_by_conf, aes(size_rt, group=conf_group, color=conf_group)) +geom_density()
# This graph would be good to present.
man_crowns_by_conf %>%
  group_by(conf_group) %>%
  summarise(mean(size_rt), mean(size_rt), median(size_rt), median(cost), sd(cost), sd(size_rt))



ggplot(man_crowns@data, aes(x=trHghts_mx, y=cost)) +
  geom_point() +
  geom_smooth(method="lm")
summary(lm(cost~trHghts_mx, data=man_crowns@data))

sum(man_crowns$shpbnd>=shpbnd_thresh)
hist(man_crowns$cost[man_crowns$shpbnd>=shpbnd_thresh]); median(man_crowns$cost[man_crowns$shpbnd>=shpbnd_thresh])
hist(man_crowns$size_rt[man_crowns$shpbnd>=shpbnd_thresh]); median(man_crowns$size_rt[man_crowns$shpbnd>=shpbnd_thresh])


man_crowns_centroid<-gCentroid(man_crowns, byid = TRUE, id = as.character(man_crowns$id))
man_crowns_coords<-coordinates(man_crowns_centroid)
xy_cost<-data.frame(x=man_crowns_coords[,'x'],
                    y=man_crowns_coords[,'y'],
                    z=man_crowns$cost)

fm.cost1 <-mgcv::gam(z ~ s(x, y, k=30),
                    data= xy_cost)
summary(fm.cost1)


man_conf_crowns<-man_crowns[man_crowns$shpbnd>=0.8,]

ggplot(man_conf_crowns@data, aes(x=size_rt)) + geom_density()

size_ration_boot<-lapply(1:1000, function(x) quantile(sample(man_conf_crowns$size_rt, size = 30, replace = TRUE), probs = c(0.3, 0.5, 0.7)))
size_ration_boot<-do.call(rbind, size_ration_boot)
apply(size_ration_boot, 2, mean)

