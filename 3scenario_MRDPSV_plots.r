# Plot boxplots of gamma (phylogenetic tree shape statistic) across 3 scenarios
# (K gradient, disturbance gradient, or no limit) and 2 regions of origin.

Allen = 0;

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

which.sims = c(3465:4064); length(which.sims);

#rootclade.stats = compile.firstlines(sim_dir,"SENC_Stats_sim")
simstats = read.csv(paste(analysis_dir,"/rootclade_stats_sims2925-4064.csv",sep=""), header=T)
sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=""),header=T)

simstats2 = merge(simstats, sim.matrix[,c(1,3:10,15,16)],by.x='sim',by.y='sim.id',all.x=T)
simstats3 = subset(simstats2, sim %in% which.sims)

simstats3$scenario = paste(simstats3$carry.cap,simstats3$energy.gradient)
simstats3$reg.of.origin = as.character(simstats3$reg.of.origin)

tropical.shade = rgb(255, 0, 0, alpha=50, maxColorValue=255)
temperate.shade = rgb(0, 0, 255, alpha=50, maxColorValue=255)

pdf(paste(analysis_dir,'/summaryplots/3scenarios_MRDPSV_',Sys.Date(),'.pdf',sep=''),height=6,width=8)
par(mfrow=c(2,1), mar=c(1,4,1,1), oma=c(2,1,1,1), las=1)

boxplot(simstats3$r.MRD.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="MRD-Richness correlation", xaxt="n",
        col = 'white',lwd=0.1,cex.axis=0.8)
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$r.MRD.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="", xaxt="n", yaxt="n",
        col = 'white',add=T,lwd=0.1,cex.lab=1,cex.axis=.8)

mtext("Temperate Origin",side=3,adj=0.2,line=0.5,cex=1)
mtext("Tropical Origin",side=3,adj = 0.8, line=0.5,cex=1)

boxplot(simstats3$r.PSV.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="PSV-Richness correlation", xaxt="n",
        col = 'white',lwd=0.1,cex.axis=0.8)
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$r.PSV.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="", xaxt="n", yaxt="n",
        col = 'white',add=T,lwd=0.1,cex.lab=1,cex.axis=0.8)

#boxplot(simstats3$r.env.MRD ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.MRD", xaxt="n",
#        col = c(rep('white',3),rep('gray50',3)))
#  rect(0, -999, 3.5, 999,col=temperate.shade)
#  rect(3.5, -999, 7, 999,col=tropical.shade)
#  boxplot(simstats3$r.env.MRD ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.MRD", xaxt="n",
#        col = c(rep('white',3),rep('gray50',3)),add=T)

#boxplot(simstats3$r.env.PSV ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.PSV", xaxt="n",
#        col = c(rep('white',3),rep('gray50',3)))
#  rect(0, -999, 3.5, 999,col=temperate.shade)
#  rect(3.5, -999, 7, 999,col=tropical.shade)
#  boxplot(simstats3$r.env.PSV ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.PSV", xaxt="n",
#        col = c(rep('white',3),rep('gray50',3)),add=T)

axis(1,rep(c('Time','Disturbance',''),2),at=1:6,cex.axis=1)
mtext(rep('Energetic\nGradient',2), 1,at=c(3,6), line=1.5)

dev.off()
