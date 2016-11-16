require(quantreg)
require(dplyr)

b<-read.csv('/Users/Tom/Dropbox/Tree segmentation/Global tree allometries/GlobalAllometricDatabase/Data-Table 1.csv')
b<-subset(b, Biogeographic_zone=="Indo-Malaya" & Biome=='Tropical forests' & Functional_type=='Angio')

b$radius<-b$CD/2
colnames(b)[colnames(b)=="H"]<-"Hmax"

########################### working with 99th percentile#################

#par( mai = c(1.2,1.2,.1,.1))
#with(b, plot(log(Hmax), log(radius), typ = "n",las = 1, tcl = -.3, xlab = "log(Tree height)", ylab = "log(Crown radius)")) 
#with(b[b$Soil == "alluvial",], points(log(Hmax), log(radius), pch =16, col ="orange")) 
#with(b[b$Soil == "kerangas",], points(log(Hmax), log(radius), pch =16, col ="slateblue")) 
#with(b[b$Soil == "sandstone",], points(log(Hmax), log(radius), pch =16, col ="grey30")) 


tau=.99
htod_lut<-function(tau)
{
  ah = rq(log(radius) ~ log(Hmax),  data = b[b$Hmax>20,], tau = tau )
  aa <- coef(ah)[1]
  bb <- coef(ah)[2]

  params<-data.frame(tau=tau, aa=aa, bb=bb)
  rownames(params)<-NULL
  return(params)
}

lut<-lapply((1:100)/100, htod_lut)
lut<-do.call(rbind, lut)

col_ramp<-colorRampPalette(c("blue", "red"))( 100 )
par(mfrow=c(1,1), mar=c(4,4,2,2))
plot(log(b$Hmax), log(b$radius),pch=19, xlab = "log(Tree height)", ylab = "log(Crown radius)",cex.axis=1.5,cex.lab=1.5)
for(i in 1:100)
  abline(lut$aa[i], lut$bb[i], col=col_ramp[i])


write.csv(lut, 'htor_lut_gta.csv')
