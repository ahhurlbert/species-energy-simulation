# Function for making box plots of the env-richness and time-richness correlation coefs across scenarios

Allen = 0;

if (Allen == 1) {
  repo_dir = 'C:/Documents and Settings/Hurlbert/species-energy-simulation'
  sim_dir = 'C:SENCoutput'
}

if (Allen ==0) {
  
  repo_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation'
  sim_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204'
  
}

sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=''),header=T)

sim.df = data.frame(K.sims.trop = 3325:3334, K.sims.temp=3345:3354, D.sims.trop = 3445:3454, D.sims.temp=3455:3464, T.sims.trop=3365:3374, T.sims.temp=3385:3394)

stats.output.tmp = c()

for (i in 1:6) {
  for (j in 1:10) {
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

write.csv(stats.output.tmp, paste(sim_dir,'/SENC_Stats_4_corr_boxplots.csv',sep=''), row.names=F)

stats.output.tmp$reg.of.origin = as.character(stats.output.tmp$reg.of.origin);
stats.output.tmp = as.data.frame(stats.output.tmp);

min.num.data.pts = 10
min.num.spp.per.clade = 30
min.num.regions = 5

tropical.shade = rgb(255, 0, 0, alpha=50, maxColorValue=255)
temperate.shade = rgb(0, 0, 255, alpha=50, maxColorValue=255)

par(mfcol=c(2,1), mar = c(4,4,1,1), oma = c(2,1,1,1))
  

boxplot(stats.output.tmp$r.env.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="",
        col = c(rep('white',3),rep('gray50',3)))
rect(0, -999, 3.5, 999,col=temperate.shade)
rect(3.5, -999, 7, 999,col=tropical.shade)
boxplot(stats.output.tmp$r.env.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="Environment-Richness Correlation", xaxt="n",
        col = c(rep('white',3),rep('gray50',3)),add=T)

mtext("Bin 10 Origin (Temperate)",side=3,adj=0,line=0.5,cex=0.375)
mtext("Bin 1 Origin (Tropical)",side=3,adj = 1, line=0.5,cex=0.375)

boxplot(stats.output.tmp$r.time.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="",
        col = c(rep('white',3),rep('gray50',3)))
rect(0, -999, 3.5, 999,col=temperate.shade)
rect(3.5, -999, 7, 999,col=tropical.shade)
boxplot(stats.output.tmp$r.time.rich ~ stats.output.tmp$scenario + stats.output.tmp$reg.of.origin, ylab="Time-Richness Correlation", 
        col = c(rep('white',3),rep('gray50',3)),add=T)




