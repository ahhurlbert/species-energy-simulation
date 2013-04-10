# Plot phylogenetic gamma statistic as a function of clade richness for the 3 different scenarios

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


Dtrop = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.trop.csv',sep=''), header=T)
Ttrop = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.csv',sep=''), header=T)
Ktrop = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)

Dtemp = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.temp.csv',sep=''), header=T)
Ttemp = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.temp.csv',sep=''), header=T)
Ktemp = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)

pdf(paste(analysis_dir,'/summaryplots/3scenarios_gamma_claderich.pdf',sep=''),height=5, width=10)
par(mfrow=c(1,2), mar = c(4,4,2,1))
cexpt = .8
plot(Ttrop$clade.rich, Ttrop$gamma.stat, ylim=range(c(Dtrop$gamma.stat, Ttrop$gamma.stat, Ktrop$gamma.stat)),
     xlim = c(0,800),pch=16, col='olivedrab3',xlab="Clade Richness",ylab="Gamma", las=1,
     main='Tropical Origin', cex = cexpt)

points(Dtrop$clade.rich, Dtrop$gamma.stat, pch = 15, col = 'salmon2', cex = cexpt)
points(Ktrop$clade.rich, Ktrop$gamma.stat, pch=17, col='skyblue', cex = cexpt)

legend('bottomleft',c('Disturbance','Time','Energy Gradient'),pch=15:17, col=c('salmon2','olivedrab3','skyblue'))

plot(Ttemp$clade.rich, Ttemp$gamma.stat, ylim=range(c(Dtrop$gamma.stat, Ttrop$gamma.stat, Ktrop$gamma.stat)),
     xlim = c(0,800),pch=16, col='olivedrab3',xlab="Clade Richness",ylab="Gamma", las=1, 
     main='Temperate Origin', cex = cexpt)

points(Dtemp$clade.rich, Dtemp$gamma.stat, pch = 15, col = 'salmon2', cex = cexpt)
points(Ktemp$clade.rich, Ktemp$gamma.stat, pch=17, col='skyblue', cex = cexpt)

dev.off()
