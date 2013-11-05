# Plot boxplots of gamma (phylogenetic tree shape statistic) across 3 scenarios
# (K gradient, disturbance gradient, or no limit) and 2 regions of origin.

Allen = 0;

if (Allen ==1) {
  
  sim_dir = "C:/SENCoutput/senc_reps_analysis"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses"
  repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"

}

if (Allen == 0) {
  
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  repo_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation"  
  
}

library(ape);

# Sebastes data
# MRD-PSV-Richness analyses
phy = read.tree(paste(repo_dir,'/Sebastes_tree_Ingram2011PRSB.phy',sep=''))
sebastes = read.csv(paste(repo_dir,'/sebastes_data_for_allen.csv',sep=''),header=T)
#Drop non-NEP species (with no latitude data)
nonNEPsp = as.character(sebastes[is.na(sebastes$min_latitude), 'X'])
NEPphy = drop.tip(phy,nonNEPsp)

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

# For Energy Gradient temperate origin,
#   MRD-S correlation predicted to be positive 
#   PSV-S correlation predicted to be negative
cor(output2)

#restricting analysis to north of Point Conception
output3 = output2[output2$lat >= 34,]
cor(output3)


# Simulation data analysis
which.sims = c(3465:4064); length(which.sims);

#rootclade.stats = compile.firstlines(sim_dir,"SENC_Stats_sim")
simstats = read.csv(paste(analysis_dir,"/rootclade_stats_sims2925-4064.csv",sep=""), header=T);
simstats$timeslice = 'no';
slicestats = read.csv(paste(analysis_dir,"/rootclade_slice_stats_sims3665-3861.csv",sep=""), header=T);
slicestats$timeslice = 'yes';
simstats = rbind(simstats,slicestats);

sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=""),header=T)

simstats2 = merge(simstats, sim.matrix[,c(1,3:10,15,16)],by.x='sim',by.y='sim.id',all.x=T)
simstats3 = subset(simstats2, sim %in% which.sims)

simstats3$scenario = NA
simstats3$scenario[simstats3$carry.cap=="on" & simstats3$energy.gradient=="on" & simstats3$timeslice=="no"] = "1 energy gradient"
simstats3$scenario[simstats3$carry.cap=="off" & simstats3$energy.gradient=="off" & simstats3$timeslice=="no"] = "2 no zero sum"
simstats3$scenario[simstats3$carry.cap=="on" & simstats3$energy.gradient=="on" & simstats3$timeslice=="yes"] = "3 pre-equilibrium"
simstats3$reg.of.origin = as.character(simstats3$reg.of.origin)

tropical.shade = 'red'
temperate.shade = 'blue'

boxplot.cols = c(rep(temperate.shade,3),rep(tropical.shade,3))

pdf(paste(analysis_dir,'/summaryplots/3scenarios_MRDPSV_',Sys.Date(),'.pdf',sep=''),height=6,width=8)
par(mfrow=c(2,1), mar=c(1,5,1,1), oma=c(5,1,1,1), las=1)

boxplot(simstats3$r.MRD.rich ~ simstats3$scenario + simstats3$reg.of.origin, 
        ylab=expression(italic(r)[MRD-richness]), xaxt="n",
        col = boxplot.cols, border = boxplot.cols, lwd=1,cex.axis=0.8, cex.lab = 1.5)
rect(0, 0, 7, 1.1, col = rgb(.6,.6,.6,.1), border = NA)
abline(h = cor(output2$MRD, output2$S), lwd = 2) #entire gradient
abline(h = cor(output3$MRD, output3$S), lwd = 2, lty = 'dashed') #north of 34N

mtext("(a)", 2, at = 1.2, las = 1, line = 3, cex = 2)

boxplot(simstats3$r.PSV.rich ~ simstats3$scenario + simstats3$reg.of.origin, 
        ylab=expression(italic(r)[PSV-richness]), xaxt="n",
        col = boxplot.cols, border = boxplot.cols,lwd=1,cex.axis=0.8, cex.lab=1.5)
rect(0, 0, 7, 1.1, col = rgb(.5,.5,.5,.1), border = NA)
abline(h = cor(output2$PSV, output2$S), lwd = 2) #entire gradient
abline(h = cor(output3$PSV, output3$S), lwd = 2, lty = 'dashed') #north of 34N

mtext("(b)", 2, at = 1.2, line = 3, las = 1, cex = 2)
axis(1, at = 1:6, labels = F)
mtext(rep(c('Energy','No zero','Pre-'),2), 1, at = 1:6, line = 1)
mtext(rep(c('gradient','sum','equilibrium'),2), 1, at = 1:6, line = 2)
mtext("Temperate Origin", side = 1, adj = 0.3, cex = 1.25, line = 3, outer = T)
mtext("Tropical Origin", side = 1, adj = 0.8, cex = 1.25, line = 3, outer = T)

#legend("bottom",c("Entire gradient","North of 34N"), lty = c('solid','dashed'))
dev.off()
