# Plot boxplots of gamma (phylogenetic tree shape statistic) across 3 scenarios
# (K gradient, disturbance gradient, or no limit) and 2 regions of origin.
sim_dir = "C:/SENCoutput"
analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"

#rootclade.stats = compile.firstlines(sim_dir,"SENC_Stats_sim")
rootclade.stats = read.csv(paste(analysis_dir,"/rootclade_stats_sims1-3464.csv",sep=""), header=T)
sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=""),header=T)

# Scenarios
dist.tr = subset(sim.matrix, max.K==22000 & reg.of.origin=='tropical')
dist.te = subset(sim.matrix, max.K==22000 & reg.of.origin=='temperate')
Kgrad.tr = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='tropical' & carry.cap=='on' & energy.gradient=='on' & max.richness==10000)
Kgrad.te = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='temperate' & carry.cap=='on' & energy.gradient=='on' & max.richness==10000)
NoLim.tr = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='tropical' & carry.cap=='off' & energy.gradient=='off' & max.richness==10000)
NoLim.te = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='temperate' & carry.cap=='off' & energy.gradient=='off' & max.richness==10000)

sub.sim = rbind(Kgrad.tr, Kgrad.te, dist.tr, dist.te, NoLim.tr, NoLim.te)

sub.gam = merge(sub.sim[,1:17], rootclade.stats[,c('sim','gamma.stat')], by.x='sim.id', by.y='sim', all.x=T)
sub.gam$scenario = paste(sub.gam$carry.cap,sub.gam$energy.gradient)
sub.gam$reg.of.origin = as.character(sub.gam$reg.of.origin)

pdf(paste(analysis_dir,'/summaryplots/3scenarios_gamma.pdf',sep=''),height=6,width=6)
par(mfrow=c(1,1), mar=c(3,4,1,1), oma=c(2,1,1,1), las=1)
boxplot(sub.gam$gamma.stat ~ sub.gam$reg.of.origin + sub.gam$scenario, ylab="Gamma", xaxt="n",
        col = rep(c('gray50','white'),3))
axis(1,c('Time', 'Disturbance','K gradient'),at=c(1.5,3.5,5.5))
legend("bottomleft",c('temperate origin','tropical origin'), pch=22, pt.bg=c('gray50','white') , pt.cex=2)
mtext("Simulation Scenario",1,adj = .55, outer=T)
dev.off()
