Load in the manual tree crowns from the lidar and UAV calculate the max heights

```{r manual crowns, include = FALSE, }
man_lid<-readOGR("../../data/shape/Manual trees LiDAR UTM.shp")
man_lid$Area<-as.numeric(as.character(man_lid$Area))
man_lid$R<-sqrt(man_lid$Area/pi) # Calculate a pseudoradius
man_lid$lid_CH_mean<-raster::extract(lid_chm, man_lid, fun=mean)
man_lid$lid_CH_max<-raster::extract(lid_chm, man_lid, fun=max)

man_uav<-readOGR("../../data/shape/Manual trees UTM.shp")
man_uav$Area<-as.numeric(as.character(man_uav$Area))
man_uav$R<-sqrt(man_uav$Area/pi) # Calculate a pseudoradius
man_uav$uav_CH_mean<-raster::extract(uav_chm, man_uav, fun=mean)
man_uav$uav_CH_max<-raster::extract(uav_chm, man_uav, fun=max)

man_CH_max_diff<-man_lid$lid_CH_max- man_uav$uav_CH_max
hist(man_CH_max_diff)
mean(man_CH_max_diff)
sd(man_CH_max_diff)

```

A comparison of LiDAR and UAV tree heights and some models to compare possible corrections.

```{r}
# Create a data frame of tree heights from manually segemented trees:
heights<-data.frame(lid = man_lid$lid_CH_max, uav = man_uav$uav_CH_max)

# K-fold cross-validation
# Fit lm model using 10-fold CV: model
intrain<-createDataPartition(y = heights$uav, times = 1, p = 0.6, list = FALSE)
training<-heights[intrain[,1],]
testing<-heights[-intrain[,1],]
fm <- train(
  lid ~ uav, training,
  method = "lm",
  trControl = trainControl(
    method = "cv", number = 10,
    verboseIter = TRUE)
  #tuneGrid  = expand.grid(intercept = FALSE)
)
summary(fm)

# Fit lm model using 10-fold CV: model forcing intercept through 0
fm.no_int <- train(
  lid ~ uav, training,
  method = "lm",
  trControl = trainControl(
    method = "cv", number = 10,
    verboseIter = TRUE),
  tuneGrid  = expand.grid(intercept = FALSE)
)
summary(fm.no_int)

# Add the prediction to the test data:
testing$pred_fm<-predict(fm, newdata = testing)
testing$pred_fm.no_int<-predict(fm.no_int, newdata = testing)


ggplot(heights, aes(x=uav, y = lid)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_abline(intercept = 0, slope = coef(fm.no_int$finalModel), color = "green") +
  geom_abline(intercept = coef(fm$finalModel)[1], slope = coef(fm$finalModel)[2], color = "blue") +
  geom_abline(intercept = 4.14, slope = 1, color = "red") +

  ggplot(heights, aes(x=lid, y = uav)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = "black") +
  geom_smooth(method = "lm", formula = y~(x)) +
  xlim(0,50)


```





Calculate the bootstrapped mean difference in LiDAR and UAV measured tree heights for a sample of n trees. This is used to calculate the coefficient of variance in the estimate of the mean difference for the sample size.

```{r bootstrap difference estimate, include = FALSE, }

ntree<-seq(10, 500, by = 10)

boot_ntree<-lapply(ntree, function(i) {
  man_CH_max_boot <- sapply(1:2000, function(j) {
    mysample <- sample(man_CH_max_diff, size = i, replace = TRUE)
    return(mean(mysample))

  })
  cov1<-(sd(man_CH_max_boot))/mean(man_CH_max_boot)
  cov2<-(sd(man_CH_max_boot)*2)/mean(man_CH_max_boot)
  return(c(cov1, cov2))
})

boot_ntree<-do.call(rbind, boot_ntree)
boot_ntree<-rbind(data.frame(n = ntree, type = 1, cov = boot_ntree[,1]),
                  data.frame(n = ntree, type = 2, cov = boot_ntree[,2])
)

```

```{r}
gg_boot_ntree<-ggplot(boot_ntree, aes(n, cov, color = as.factor(type))) +
  geom_line()
print(gg_boot_ntree)
```


This is interesting because it suggests that if we sample 100 trees, we will be within 15% of the true value 68% of the time and 30% of the true value 95% of the time.


```{r difference sample of 100}
man_CH_max_boot <- sapply(1:2000, function(i) {
  mysample <- sample(man_CH_max_diff, size = 30, replace = TRUE)
  return(mean(mysample))
})

hist(man_CH_max_boot)
(quants<-quantile(man_CH_max_boot, probs = c(0.05, 0.25, 0.5, 0.75, 0.95)))
```
