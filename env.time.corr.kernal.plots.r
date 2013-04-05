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

stats.output.tmp = read.csv(paste(sim_dir,'/SENC_Stats_4_corr_boxplots.csv',sep=''));

stats.output.tmp$reg.of.origin = as.character(stats.output.tmp$reg.of.origin);
stats.output.tmp = as.data.frame(stats.output.tmp);

min.num.data.pts = 10
min.num.spp.per.clade = 30
min.num.regions = 5

stats.output.tmp = stats.output.tmp[which(stats.output.tmp$n.regions >= min.num.regions),];
stats.output.tmp = stats.output.tmp[which(stats.output.tmp$extant.richness >= min.num.spp.per.clade),];
stats.output.tmp = stats.output.tmp[which(is.na(stats.output.tmp$r.time.rich)==F),];
stats.output.tmp = stats.output.tmp[which(is.na(stats.output.tmp$r.env.rich)==F),];

head(stats.output.tmp); dim(stats.output.tmp);
range(stats.output.tmp$n.regions);
range(stats.output.tmp$extant.richness);

tropical.shade = rgb(255, 0, 0, alpha=50, maxColorValue=255)
temperate.shade = rgb(0, 0, 255, alpha=50, maxColorValue=255)

band = 0.05;

pdf(paste(fig_dir,"/SENC_r.E.R_r.T.R._kernals.pdf",sep=""),height=7,width=7)

  par(mfrow=c(2,2), mar = c(4,4,1,1), oma = c(2,1,1,1))
  
  plot(density(stats.output.tmp$r.env.rich[which(is.element(stats.output.tmp$sim,sim.df[,'K.sims.temp'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,ylim=c(0,7),xlim=c(-1,1),main="",xlab="r (Env-Richness)",col=2)
  points(density(stats.output.tmp$r.env.rich[which(is.element(stats.output.tmp$sim,sim.df[,'D.sims.temp'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,col=4)
  points(density(stats.output.tmp$r.env.rich[which(is.element(stats.output.tmp$sim,sim.df[,'T.sims.temp'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,lty=2,col=3)
  text(-1.1,6.9,labels="Temperate Origin",pos=4)

  plot(density(stats.output.tmp$r.env.rich[which(is.element(stats.output.tmp$sim,sim.df[,'K.sims.trop'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,ylim=c(0,7),xlim=c(-1,1),main="",xlab="r (Env-Richness)",col=2)
  points(density(stats.output.tmp$r.env.rich[which(is.element(stats.output.tmp$sim,sim.df[,'D.sims.trop'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,col=4)
  points(density(stats.output.tmp$r.env.rich[which(is.element(stats.output.tmp$sim,sim.df[,'T.sims.trop'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,lty=2,col=3)
  text(-1.1,6.9,labels="Tropical Origin",pos=4)
  
  plot(density(stats.output.tmp$r.time.rich[which(is.element(stats.output.tmp$sim,sim.df[,'K.sims.temp'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,ylim=c(0,7),xlim=c(-1,1),main="",xlab="r (Time-Richness)",col=2)
  points(density(stats.output.tmp$r.time.rich[which(is.element(stats.output.tmp$sim,sim.df[,'D.sims.temp'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,col=4)
  points(density(stats.output.tmp$r.time.rich[which(is.element(stats.output.tmp$sim,sim.df[,'T.sims.temp'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,lty=2,col=3)
  text(-1.1,6.9,labels="Temperate Origin",pos=4)
  legend(-1,6,legend=c("Time",'Disturbance','Energy Gradient'),lty=c(2,1,1),col=c(3,2,4),lwd=2)

  plot(density(stats.output.tmp$r.time.rich[which(is.element(stats.output.tmp$sim,sim.df[,'K.sims.trop'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,ylim=c(0,7),xlim=c(-1,1),main="",xlab="r (Time-Richness)",col=2)
  points(density(stats.output.tmp$r.time.rich[which(is.element(stats.output.tmp$sim,sim.df[,'D.sims.trop'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,col=4)
  points(density(stats.output.tmp$r.time.rich[which(is.element(stats.output.tmp$sim,sim.df[,'T.sims.trop'])==T)],bw=band,from=-1,to=1),typ="l",lwd=2,lty=2,col=3)
  text(-1.1,6.9,labels="Tropical Origin",pos=4)

dev.off();
