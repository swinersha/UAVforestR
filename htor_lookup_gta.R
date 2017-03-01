require(quantreg)
require(dplyr)

source("R/htor_lut.R") # The funciton to run the quantile regression

b<-read.csv('/Users/Tom/Dropbox/Tree segmentation/Global tree allometries/GlobalAllometricDatabase/Data-Table 1.csv')
b<-subset(b, Biogeographic_zone=="Indo-Malaya" & Biome=='Tropical forests' & Functional_type=='Angio')

b$radius<-b$CD/2
colnames(b)[colnames(b)=="H"]<-"Hmax"

# Run the quantile regressions:
lut<-htod_lut(height=b$Hmax, radius=b$radius)

col_ramp<-colorRampPalette(c("blue", "red"))( 100 )
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(log(b$Hmax), log(b$radius),pch=19, xlab = "log(Tree height)", ylab = "log(Crown radius)",cex.axis=1.5,cex.lab=1.5)
for(i in 1:100)
  abline(lut$aa[i], lut$bb[i], col=col_ramp[i])


write.csv(lut, 'htor_lut_gta.csv')
