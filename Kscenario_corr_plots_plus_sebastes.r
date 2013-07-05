# Figure 4: Plot latitude-richness correlation versus (a) time of clade origin and
# (b) clade richness for all simulation subclades with at least 30 species and
# spanning at least 5 regions. (c) Plot same thing for empirical Sebastes data
# for the entire northeastern Pacific gradient (23-66N) as well as the gradient
# north of Point Conception (34-66N).

# libraries
require(ape)
require(caper)

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

#Simulation data
Ktrop.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)
Ktemp.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)

Ktrop = subset(Ktrop.all, clade.richness >= 30 & n.regions >=5)
Ktemp = subset(Ktemp.all, clade.richness >= 30 & n.regions >=5)

#Sebastes data
sebastes = read.csv(paste(repo_dir,'/sebastes_data_for_allen.csv',sep=''),header=T)
phy = read.tree(paste(repo_dir,'/Sebastes_tree_Ingram2011PRSB.phy',sep=''))
#Drop non-NEP species (with no latitude data)
nonNEPsp = as.character(sebastes[is.na(sebastes$min_latitude), 'X'])
NEPphy = drop.tip(phy,nonNEPsp)
lat = min(sebastes$min_latitude, na.rm=T):max(sebastes$max_latitude, na.rm=T)

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



#Plot
pdf(paste(analysis_dir,'/summaryplots/latcorr_subclades_vs_cladeAgeRich_',Sys.Date(),'.pdf',sep=''), height=15, width = 7)
par(mfrow=c(3,1), mar = c(6, 3, 2, 1), oma=c(1, 10, 1, 1), mgp = c(5, 1.5, 0))
cexpts = 2
cexpts.seb = 3
cexaxis = 2
cexlab = 3
cexleg = 2
cexabc = 2.5
pch.temp = 18

#Vs clade origin time (negative the r.env.rich is equal to the richness-latitude correlation)
plot(log10(Ktrop$clade.origin.time), -Ktrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),
     xlab = expression(paste(plain(log)[10]," Time of Clade Origin")),
     ylab = "", cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain, las=1)
points(log10(Ktemp$clade.origin.time), -Ktemp$r.env.rich, col = 'blue', pch = pch.temp, cex = cexpts)
abline(h = 0, lty = 'dashed')
legend('topleft',c('temperate origin','tropical origin'), pch = c(pch.temp, 16), col = c('blue','red'), cex = cexleg)
mtext("(a)", 2, at = 1, cex = cexabc, outer = T, las = 1, line = 5)

#Vs clade richness
plot(log10(Ktrop$clade.richness), -Ktrop$r.env.rich, pch = 16, col = 'red', ylim = c(-1,1),
     xlab = expression(paste(plain(log)[10]," Clade Richness")), ylab = "",
     cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain, las = 1)
points(log10(Ktemp$clade.richness), -Ktemp$r.env.rich, col = 'blue', pch = pch.temp, cex = cexpts)
abline(h = 0,lty = 'dashed')
mtext("Latitude-richness correlation", 2, outer=T, cex = cexaxis, line = 3)
mtext("(b)", 2, outer=T, at = .67, cex = cexabc, las = 1, line = 5)

#extra tick marks showing % of max richness
#pcts = c(.9,.3,.1,.03,.01)
#par(mgp=c(-3,0,0))
#axis(1,at=log10(pcts*max(Ktrop$clade.richness)), labels=F,tck= .01)
#text(log10(pcts*max(Ktrop$clade.richness)), rep(-1,length(pcts)), paste(pcts*100,"%",sep=""))

#Sebastes vs clade richness
plot(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich, ylim = c(-1,1), las=1, cex.axis = cexaxis, cex.lab = cexlab,
     xlab = expression(paste(plain(log)[10]," Clade Richness")), ylab = 'Latitude-richness correlation', pch = 16, cex = cexpts.seb)
points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.gte34, pch = 1, cex = cexpts.seb)
#points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.lt34, pch=2, cex = cexpts)
abline(h=0,lty='dashed')
legend("topright", c('Entire gradient','North of 34N'), pch = c(16, 1), cex = cexleg)
mtext("(c)", 2, at = 0.33, outer = T, cex = cexabc, las = 1, line = 5)

dev.off()


