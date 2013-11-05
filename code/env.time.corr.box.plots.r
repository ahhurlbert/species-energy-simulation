# Function for making box plots of the env-richness and time-richness correlation coefs across scenarios

Allen = 0;

if (Allen == 1) {
  repo_dir = 'C:/Documents and Settings/Hurlbert/species-energy-simulation'
  sim_dir = 'C:SENCoutput'
}

if (Allen ==0) {
  
  repo_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation'
  sim_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204'
  fig_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204/summaryplots'
  
}

sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=''),header=T)

sim.df = data.frame(K.sims.trop = 3665:3764, K.sims.temp = 3765:3864, D.sims.trop = 3865:3964, D.sims.temp = 3965:4064, T.sims.trop=3465:3564, T.sims.temp=3565:3664)

stats.output.tmp = c()

for (i in 1:6) {
  for (j in 1:100) {
    tmp = read.csv(paste(sim_dir,'/SENC_stats_sim',sim.df[j,i],'.csv',sep=''), header=T)
    tmp$reg.of.origin = sim.matrix$reg.of.origin[sim.matrix$sim.id == sim.df[j,i]]
    if (sim.matrix$carry.cap[sim.matrix$sim.id == sim.df[j,i]] == 'off') {scenario = 'Time'};
    if (sim.matrix$carry.cap[sim.matrix$sim.id == sim.df[j,i]] == 'on' & sim.matrix$energy.gradient[sim.matrix$sim.id == sim.df[j,i]] == 'off') {scenario = 'Disturbance'}
    if (sim.matrix$carry.cap[sim.matrix$sim.id == sim.df[j,i]] == 'on' & sim.matrix$energy.gradient[sim.matrix$sim.id == sim.df[j,i]] == 'on') {scenario = 'Energy Gradient'}
    tmp$scenario = scenario
    stats.output.tmp = rbind(stats.output.tmp, tmp)
    rm('scenario','tmp')
  }
  
}

head(stats.output.tmp); dim(stats.output.tmp);

write.csv(stats.output.tmp, paste(sim_dir,'/SENC_Stats_4_corr_boxplots.csv',sep=''), row.names=F,quote=F)

stats.output.tmp$reg.of.origin = as.character(stats.output.tmp$reg.of.origin);
stats.output.tmp = as.data.frame(stats.output.tmp);

min.num.data.pts = 10
min.num.spp.per.clade = 30
min.num.regions = 5

stats.output.tmp = stats.output.tmp[which(stats.output.tmp$n.regions >= min.num.regions),];
stats.output.tmp = stats.output.tmp[which(stats.output.tmp$extant.richness >= min.num.spp.per.clade),];
head(stats.output.tmp); dim(stats.output.tmp);
range(stats.output.tmp$n.regions);
range(stats.output.tmp$extant.richness);

tropical.shade = rgb(255, 0, 0, alpha=50, maxColorValue=255)
temperate.shade = rgb(0, 0, 255, alpha=50, maxColorValue=255)

pdf(paste(fig_dir,"/SENC_r.E.R_r.T.R._boxplots.pdf",sep=""),height=7,width=7)

  par(mfcol=c(2,1), mar = c(4,4,1,1), oma = c(2,1,1,1))
  

  boxplot(stats.output.tmp$r.env.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="",
        col = c(rep('white',3),rep('gray50',3)),xaxt="n")
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(stats.output.tmp$r.env.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="Environment-Richness Correlation", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T,xaxt="n")

  mtext("Bin 10 Origin (Temperate)",side=3,adj=0,line=0.5,cex=1)
  mtext("Bin 1 Origin (Tropical)",side=3,adj = 1, line=0.5,cex=1)

  boxplot(stats.output.tmp$r.time.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="",
        col = c(rep('white',3),rep('gray50',3)),xaxt="n")
  rect(0, -999, 3.5, 999,col=temperate.shade)
  rect(3.5, -999, 7, 999,col=tropical.shade)
  boxplot(stats.output.tmp$r.time.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="Time-Richness Correlation", 
        col = c(rep('white',3),rep('gray50',3)),add=T,xaxt="n")

  axis(1,rep(c('D','E.G.','T'),2),at=c(1:6))

dev.off();


