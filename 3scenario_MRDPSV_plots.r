# Plot boxplots of gamma (phylogenetic tree shape statistic) across 3 scenarios
# (K gradient, disturbance gradient, or no limit) and 2 regions of origin.
sim_dir = "C:/SENCoutput"
analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"

which.sims = c(3325:3334,3345:3354,3365:3374,3385:3394,3445:3464)

#rootclade.stats = compile.firstlines(sim_dir,"SENC_Stats_sim")
simstats = read.csv(paste(analysis_dir,"/rootclade_stats_sims1-3464.csv",sep=""), header=T)
sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=""),header=T)

simstats2 = merge(simstats, sim.matrix[,c(1,3:10,15,16)],by.x='sim',by.y='sim.id',all.x=T)
simstats3 = subset(simstats2, sim %in% which.sims)

simstats3$scenario = paste(simstats3$carry.cap,simstats3$energy.gradient)
simstats3$reg.of.origin = as.character(simstats3$reg.of.origin)

pdf(paste(analysis_dir,'/summaryplots/3scenarios_MRDPSV.pdf',sep=''),height=6,width=6)
par(mfrow=c(4,1), mar=c(3,4,1,1), oma=c(2,1,1,1), las=1)
boxplot(simstats3$r.MRD.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.MRD.rich", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
boxplot(simstats3$r.PSV.rich ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.PSV.rich", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
boxplot(simstats3$r.env.MRD ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.MRD", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
boxplot(simstats3$r.env.PSV ~ simstats3$scenario + simstats3$reg.of.origin, ylab="r.env.PSV", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)))
axis(1,rep(c('Time', 'Disturbance','K gradient'),2),at=c(1:6))
mtext("Temperate origin",1,adj = .25, outer=T)
mtext("Tropical origin",1,adj = .75, outer=T)
dev.off()
