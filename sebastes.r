# Code for analyzing patterns in geographic dataset of Northeast Pacific rockfish,
# genus Sebastes, from data provided by Travis Ingram and used in his 2011 PRSB paper.

# Purpose: To estimate species richness, mean root distance, and phylogenetic species
# variability in each latitudinal bin, and to relate those values to each other.

# Author: Allen Hurlbert

setwd('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/sebastes')
require(ape)
require(plyr)
require(picante)
require(caper) #for clade.members()

sebastes = read.csv('sebastes_data_for_allen.csv',header=T)
phy = read.tree('Sebastes_tree_Ingram2011PRSB.phy')

#Drop non-NEP species (with no latitude data)
nonNEPsp = as.character(sebastes[is.na(sebastes$min_latitude), 'X'])
NEPphy = drop.tip(phy,nonNEPsp)


lat = min(sebastes$min_latitude, na.rm=T):max(sebastes$max_latitude, na.rm=T)

# Marine NPP data
npp.1dg = read.csv('pacificnppvslatitude_1dg.csv', header=T)
npp.2dg = read.csv('pacificnppvslatitude_2dg.csv', header=T)
npp.1d.ma = read.csv('pacificnppvslatitude_1dg_ma.csv', header=T) # 1deg resolution of values representing a 5 deg moving average

##############################################################################
# MRD-PSV-Richness analyses
richness = sapply(lat, function(x) nrow(subset(sebastes, min_latitude <= x & max_latitude >= x)))

phylo.bl1 <- compute.brlen(NEPphy, 1)
all.dist <- dist.nodes(phylo.bl1)
root.dist <- all.dist[length(NEPphy$tip.label)+1, 1:length(NEPphy$tip.label)]
tips.to.root <- data.frame(spp.name=NEPphy$tip.label,root.dist)

output = c()
for (i in lat) {
  species = subset(sebastes, min_latitude <= i & max_latitude >= i, select='X')
  
  #MRD
  MRD.ini <- merge(species, tips.to.root, by.x="X", by.y="spp.name",sort = FALSE)
  MRD <- mean(MRD.ini$root.dist)
  
  #PSV
  Vmatrix = vcv(NEPphy, corr=F)
  psvs = matrix(NA, ncol=2)
  
  index = row.names(Vmatrix) %in% species$X
  v.matrix = Vmatrix[index,index]
  n = nrow(v.matrix)
  psv = (n*sum(diag(v.matrix)) - sum(v.matrix))/(sum(diag(v.matrix))*(n-1))
  
  output = rbind(output, c(i, MRD, psv))
}

output2 = data.frame(cbind(output, richness))
names(output2) = c('lat','MRD','PSV','S')

output2 = merge(output2, npp.1dg, by.x='lat', by.y='latitude', all.x=T)
output2 = merge(output2, npp.2dg, by.x='lat', by.y='latitude', all.x=T)

# For Energy Gradient temperate origin,
#   MRD-S correlation predicted to be positive 
#   PSV-S correlation predicted to be negative
cor(output2, use = "pairwise.complete.obs")

#restricting analysis to north of Point Conception
output3 = output2[output2$lat >= 34,]
cor(output3)


# Supplemental figure showing species richness, MRD, and PSV of Pacific rockfish clade
# as a function of latitude
pdf(paste('sebastes_MRD-PSV_corrs_',Sys.Date(),'.pdf',sep=''),height=6,width=9)
par(mar = c(6,12,1,9), mgp = c(4.5, 2, 0))
plot(lat,richness, type = 'l', lwd = 3, xlab = "Latitude", ylab = "Richness", cex.lab = 1.5, cex.axis = 1.5, las = 1)
rect(20,-1,34,60, col = rgb(.9,.9,.9,.4), border=NA)

par(new=T, mgp = c(4,1,0))

plot(lat, output2$MRD/max(root.dist), col='blue',xaxt="n",yaxt="n",ylab="", xlab = "", pch = 16)
axis(4, line = .5, col = 'blue', las = 1, cex.axis = 1.5, lwd = 3)
par(new=T)
plot(lat,output2$PSV, col='red',xaxt="n",yaxt="n",ylab="", xlab = "", pch = 16)
axis(4, line = 5, col = 'red', las = 1, cex.axis = 1.5, lwd = 3)

par(new=T)

plot(npp.1d.ma$latitude, npp.1d.ma$npp, type='l', col='green', xlim = c(min(lat),max(lat)),
     xaxt="n", yaxt="n", lwd = 2, xlab = "", ylab="")
axis(2, line = 7, col = 'green', cex.axis = 1.5, lwd =2)

mtext(expression(paste("NPP (mg C ", m^-2, d^-1,")", sep=" ")), 2, cex = 1.5, line = 9.5)

legend("topright",c('richness','MRD','PSV','NPP'),lty = c('solid','dotted','dotted','solid'),
       col=c('black','blue','red','green'), seg.len = 1, lwd = c(4, 8, 8, 3), cex = 1)
dev.off()



############################################################################
#Gamma plot

Allen = 1;

if (Allen ==1) {
  
  sim_dir = "C:/SENCoutput/senc_reps_analysis"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
  repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"
  
}

if (Allen == 0) {
  
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  repo_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation"  
  
}

Ttrop = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.csv',sep=''), header=T)
Ktrop = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)
Ktrop.slice = read.csv(paste(sim_dir,'/SENC_Stats_K.slice.sims.trop.csv',sep=''), header=T)

Ktemp = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)
Ktemp.slice = read.csv(paste(sim_dir,'/SENC_Stats_K.slice.sims.temp.csv',sep=''), header=T)

Tline = 'olivedrab3'
Kline = 'mediumorchid2'
Kline.slice = 'goldenrod2'

#Sebastes phylogeny has 99 species (only 66 in NEP), so pull out clades for each scenario of roughly the same size
rich=66

Ttrop99 = subset(Ttrop, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktrop99 = subset(Ktrop, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktrop.slice99 = subset(Ktrop.slice, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktemp99 = subset(Ktemp, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktemp.slice99 = subset(Ktemp.slice, clade.richness > .9*rich & clade.richness < 1.1*rich)

pdf(paste(analysis_dir, '/sebastes/sebastes_gamma.pdf', sep=''), height = 6, width = 8)
plot(density(Ttrop99$gamma.stat), col=Tline, main="", xlab="Gamma", lwd=3, xlim = c(-8,2))
points(density(Ktrop99$gamma.stat), type='l',col=Kline, lty='dotted',lwd=3)
points(density(Ktemp99$gamma.stat), type='l',col=Kline, lwd=3)
abline(v = gammaStat(NEPphy), lwd=2, lty='dashed')
legend("topleft",c('no zero-sum constraint','zero-sum w/ tropical origin', 'zero-sum w/ temperate origin', 'observed gamma'),
       col = c(Tline, Kline, Kline, 'black'), lty = c('solid', 'dotted', 'solid', 'dashed'), lwd=3)
points(density(Ktrop.slice99$gamma.stat), type= 'l', col=Kline.slice, lty='dotted', lwd=2)
points(density(Ktemp.slice99$gamma.stat), type = 'l', col=Kline.slice, lwd=2)
dev.off()


#points(density(Ktrop.slice99$gamma.stat), col=Kline.slice, lty='dotted', lwd=2)
#points(density(Ktemp.slice99$gamma.stat), col=Kline.slice, lwd=2)

#############################################################################
#Plotting depth range on phylogeny
sebastes$mean_depth = apply(sebastes[,c('min_common_depth','max_common_depth')], 1, function(x) mean(x, na.rm=T))
depth.col = colorRampPalette(colors()[c(405,431,616,619,566)]) #dark blue = deep, light blue = shallow

sebastes$mean_lat = apply(sebastes[,c('min_latitude','max_latitude')], 1, function(x) mean(x, na.rm=T))
lat.col = colorRampPalette(c('red','orange','yellow','green')) #red = low latitude, green = high

plot(phy)
#reflect continuous depth
tiplabels(pch=15,col = depth.col(100)[floor(100*sebastes$mean_depth/max(sebastes$mean_depth, na.rm=T))], adj=3.8, cex=1.25)
#reflect categorical shallow/deep
depth.threshold = 180
sebastes$shde.col[sebastes$mean_depth < depth.threshold] = colors()[405]
sebastes$shde.col[sebastes$mean_depth >= depth.threshold] = colors()[566]
tiplabels(pch=15,col = sebastes$shde.col, adj=4, cex=1.25)  
#reflect latitude
tiplabels(pch=15,col = lat.col(100)[floor(100*sebastes$mean_lat/max(sebastes$mean_lat, na.rm=T))], adj=4.2, cex=1.25)


##############################################################################
# Subclade correlations

min.num.spp = 5

lat.corr.output = c()
for (c in (NEPphy$Nnode+2):max(NEPphy$edge)) {
  
  #pull out list of species names belonging to each subclade
  sub.clade = clade.members(c, NEPphy, tip.labels=T)
  sub.populations = subset(sebastes, X %in% sub.clade);
  
  sub.richness = sapply(lat, function(x) nrow(sub.populations[sub.populations$min_latitude <= x & sub.populations$max_latitude >= x, ]))
  if(length(sub.clade) >= min.num.spp) {
    lat.corr = cor(lat[sub.richness>0], sub.richness[sub.richness>0])
    lat.corr2 = cor(lat[sub.richness>0 & lat >= 34], sub.richness[sub.richness > 0 & lat >=34])
    lat.corr3 = cor(lat[sub.richness>0 & lat < 34], sub.richness[sub.richness > 0 & lat < 34])
    lat.corr.output = rbind(lat.corr.output, c(c, length(sub.clade), lat.corr, lat.corr2, lat.corr3))
  }
}
lat.corr.output = data.frame(lat.corr.output)
names(lat.corr.output) = c('cladeID','clade.richness','r.lat.rich','r.lat.rich.gte34','r.lat.rich.lt34')  
# Correlations calculated from entire gradient (Alaska-to Baja), and for the gradient
# north of 34N (Alaska-Point Conception)

pdf(paste(analysis_dir,'/sebastes/sebastes_corrplot.pdf',sep=''),height=6,width=8)
par(mar=c(4,4,1,1), pch = 1.5)
plot(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich, ylim = c(-1,1), pch = 16,
     xlab = expression(paste(plain(log)[10]," Clade Richness")),ylab = 'Latitude-richness correlation')
points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.gte34, pch=1)
points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.lt34, pch=1, col = 'red')
abline(h=0,lty='dashed')
legend(1.6,.5, c('Entire gradient','North of 34N','South of 34N'), pch = c(16,1,1), col = c('black','black','red'))
dev.off()


# As above, but looking at NPP-richness correlations instead of latitude-richness
lat2 = 23:60 #no NPP data north of 60 deg in this case
npp = npp.1d.ma$npp[order(npp.1d.ma$lat)]
npp.corr.output = c()
for (c in (NEPphy$Nnode+2):max(NEPphy$edge)) {
  
  #pull out list of species names belonging to each subclade
  sub.clade = clade.members(c, NEPphy, tip.labels=T)
  sub.populations = subset(sebastes, X %in% sub.clade);
  
  sub.richness = sapply(lat2, function(x) nrow(sub.populations[sub.populations$min_latitude <= x & sub.populations$max_latitude >= x, ]))
  if(length(sub.clade) >= min.num.spp) {
    npp.corr = cor(npp[sub.richness>0], sub.richness[sub.richness>0])
    npp.corr2 = cor(npp[sub.richness>0 & lat2 >= 34], sub.richness[sub.richness > 0 & lat2 >=34])
    npp.corr3 = cor(npp[sub.richness>0 & lat2 < 34], sub.richness[sub.richness > 0 & lat2 < 34])
    npp.corr.output = rbind(npp.corr.output, c(c, length(sub.clade), npp.corr, npp.corr2, npp.corr3))
  }
}
npp.corr.output = data.frame(npp.corr.output)
names(npp.corr.output) = c('cladeID','clade.richness','r.npp.rich','r.npp.rich.gte34','r.npp.rich.lt34')  

pdf(paste(analysis_dir,'/sebastes_NPPcorrplot.pdf',sep=''),height=6,width=8)
par(mar=c(4,4,1,1), cex = 1.5)
plot(log10(npp.corr.output$clade.richness), npp.corr.output$r.npp.rich, ylim = c(-1,1), 
     xlab = expression(paste(plain(log)[10]," Clade Richness")),ylab = 'NPP-richness correlation', pch=16)
points(log10(npp.corr.output$clade.richness), npp.corr.output$r.npp.rich.gte34, pch=1)
points(log10(npp.corr.output$clade.richness), npp.corr.output$r.npp.rich.lt34, pch=1, col = 'red')
abline(h=0,lty='dashed')
legend('bottomright', c('Entire gradient','North of 34N','South of 34N'), pch = c(16,1,1), col = c('black','black','red'))
dev.off()

