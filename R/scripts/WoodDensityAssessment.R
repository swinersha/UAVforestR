


library("BIOMASS")
library("tidyverse")

setwd('/Users/Tom/Documents/Work/RSPB/HRF/Restoration trials/Selective thinning/Block B trial 2013/Census 2/Data')
thin<-read.csv('thin_RTA.csv')

species_strings<-stringr::str_extract_all(string = thin$nama.ilmiah, pattern = "[:alpha:]*")
genus_str<-sapply(species_strings, function(x) x[1])
species_str<-sapply(species_strings, function(x) x[3])

# Correct species binomials with BIOMASS::correctTaxo: 

#binom_correct<-correctTaxo(species = species_str[1:5], genus = genus_str[1:5]) 
#getWoodDensity(species = binom_correct$speciesCorrected, genus = binom_correct$genusCorrected) 
wd<-getWoodDensity(species = species_str, genus = genus_str) 
fm.HD<-modelHD(D = thin$dbh, H = thin$tinggi, method="log2", useWeight = TRUE) 
Hpred<-retrieveH(D = thin$dbh, model = fm.HD)

plot(y= Hpred$H, x=thin$tinggi)


# Calculated plot level AGB ----

thin$AGB<-computeAGB(D = thin$dbh, WD = wd$meanWD, H = Hpred$H)
plotAGB<-thin %>% 
  filter(sp %in% c(10)) %>%
  group_by(plot, sp) %>%
  summarise(AGB = sum(AGB), comp = comp[1], sf = sf[1]) #%>%
  #mutate(AGB.ha = AGB * sf) # total biomass scaled to per hectare values
  # ungroup() %>%
  # group_by(comp) %>%
  # summarise(AGB.ha = sum(AGB.ha)) %>%
  # ggplot(aes(AGB.ha)) + geom_histogram(, bins = 15)


# Bootstrap AGB ----
nboot<-1000
smplr<-sapply(1:nboot, function(x) sample(1:nrow(plotAGB), size = 25, replace = TRUE))
AGBboot<-apply(smplr, 2, function(x){
  sum(plotAGB$AGB[x])
})
hist(AGBboot); mean(AGBboot)

# Calculated community weighted wood density ----

weighted.mean(x = wd$meanWD, w = thin$AGB)





