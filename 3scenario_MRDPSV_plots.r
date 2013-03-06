# Plot boxplots of gamma (phylogenetic tree shape statistic) across 3 scenarios
# (K gradient, disturbance gradient, or no limit) and 2 regions of origin.

Allen = 0;

if (Allen ==1) {
  
  sim_dir = "C:/SENCoutput"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
  repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"

}

if (Allen == 0) {
  
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  repo_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation"  
  
}

which.sims = c(3325:3334,3345:3354,3365:3374,3385:3394,3445:3464)

#rootclade.stats = compile.firstlines(sim_dir,"SENC_Stats_sim")
simstats = read.csv(paste(analysis_dir,"/rootclade_stats_sims1-3464.csv",sep=""), header=T)
sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=""),header=T)

simstats2 = merge(simstats, sim.matrix[,c(1,3:10,15,16)],by.x='sim',by.y='sim.id',all.x=T)
simstats3 = subset(simstats2, sim %in% which.sims)

simstats3$scenario = paste(simstats3$carry.cap,simstats3$energy.gradient)
simstats3$reg.of.origin = as.character(simstats3$reg.of.origin)

tropical.shade = rgb(255, 0, 0, alpha=50, maxColorValue=255)
temperate.shade = rgb(0, 0, 255, alpha=50, maxColorValue=255)

pdf(paste(analysis_dir,'/summaryplots/3scenarios_MRDPSVGam.pdf',sep=''),height=6,width=2.5)
par(mfrow=c(5,1), mar=c(1,4,1,1), oma=c(2,1,1,1), las=1)

boxplot(simstats3$gamma.stat ~ simstats3$scenario + simstats3$reg.of.origin, ylab="Gamma", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$gamma.stat ~ simstats3$scenario + simstats3$reg.of.origin, ylab="Gamma", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T)

  mtext("Bin 10 Origin (Temperate)",side=3,adj=0,line=0.5,cex=0.375)
  mtext("Bin 1 Origin (Tropical)",side=3,adj = 1, line=0.5,cex=0.375)

boxplot(simstats3$r.MRD.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.MRD.rich", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$r.MRD.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.MRD.rich", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T)

boxplot(simstats3$r.PSV.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.PSV.rich", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$r.PSV.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.PSV.rich", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T)

boxplot(simstats3$r.env.MRD ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.MRD", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$r.env.MRD ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.MRD", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T)

boxplot(simstats3$r.env.PSV ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.PSV", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(simstats3$r.env.PSV ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.PSV", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T)

axis(1,rep(c('T','D','K'),2),at=c(1:6))

dev.off()
