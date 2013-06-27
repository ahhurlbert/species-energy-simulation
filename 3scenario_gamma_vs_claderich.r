# Plot phylogenetic gamma statistic as a function of clade richness for the 3 different scenarios

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


Dtrop = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.trop.csv',sep=''), header=T)
Ttrop = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.csv',sep=''), header=T)
Ktrop = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)
Ktrop.slice = read.csv(paste(sim_dir,'/SENC_Stats_K.slice.sims.trop.csv',sep=''), header=T)

Dtemp = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.temp.csv',sep=''), header=T)
Ttemp = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.temp.csv',sep=''), header=T)
Ktemp = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)
Ktemp.slice = read.csv(paste(sim_dir,'/SENC_Stats_K.slice.sims.temp.csv',sep=''), header=T)

Dcol = 'gold2'
Tcol = 'olivedrab3'
Kcol = 'mediumorchid2'
Kcol.slice = 'gold2'
Dline = 'goldenrod3'
Tline = 'darkgreen'
Kline = 'mediumorchid4'
Kline.slice = 'goldenrod3'

plot.disturbance = 0; # toggle for plotting disturbance or the time slice from K

pdf(paste(analysis_dir,'/summaryplots/3scenarios_gamma_claderich_',Sys.Date(),'.pdf',sep=''),height=5, width=10)
par(mfrow=c(1,2), mar = c(4,4,2,1))
cexpt = .8

if (plot.disturbance == 1) {

plot(Ttrop$clade.rich, Ttrop$gamma.stat, ylim=range(c(Dtrop$gamma.stat, Ttrop$gamma.stat, Ktrop$gamma.stat)),
     xlim = c(0,800),pch=16, col=Tcol,xlab="Clade Richness",ylab="Gamma", las=1,
     main='Tropical Origin', cex = cexpt)

points(Dtrop$clade.rich, Dtrop$gamma.stat, pch = 15, col = Dcol, cex = cexpt)
points(Ktrop$clade.rich, Ktrop$gamma.stat, pch=17, col=Kcol, cex = cexpt)

#splines
points(smooth.spline(Dtrop$clade.rich,Dtrop$gamma.stat,df=5),type='l',col= Dline,lwd=4)
points(smooth.spline(Ttrop$clade.rich,Ttrop$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(Ktrop$clade.rich,Ktrop$gamma.stat,df=5),type='l',col= Kline,lwd=4)

legend('bottomleft',c('Time','Disturbance','Energy Gradient'),pch=c(16,15,17), col=c(Tcol,Dcol,Kcol))

plot(Ttemp$clade.rich, Ttemp$gamma.stat, ylim=range(c(Dtrop$gamma.stat, Ttrop$gamma.stat, Ktrop$gamma.stat)),
     xlim = c(0,800),pch=16, col= Tcol, xlab="Clade Richness",ylab="Gamma", las=1, 
     main='Temperate Origin', cex = cexpt)

points(Dtemp$clade.rich, Dtemp$gamma.stat, pch = 15, col = Dcol, cex = cexpt)
points(Ktemp$clade.rich, Ktemp$gamma.stat, pch=17, col= Kcol, cex = cexpt)

#splines
points(smooth.spline(Dtemp$clade.rich,Dtemp$gamma.stat,df=5),type='l',col= Dline,lwd=4)
points(smooth.spline(Ttemp$clade.rich,Ttemp$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(Ktemp$clade.rich,Ktemp$gamma.stat,df=5),type='l',col= Kline,lwd=4)

} # end if statement for plotting disturbance

if (plot.disturbance == 0) {
  
  plot(Ttrop$clade.rich, Ttrop$gamma.stat, ylim=range(c(Ktrop.slice$gamma.stat, Ttrop$gamma.stat, Ktrop$gamma.stat)),
       xlim = c(0,800),pch=16, col=Tcol,xlab="Clade Richness",ylab="Gamma", las=1,
       main='Tropical Origin', cex = cexpt)
  
  points(Ktrop$clade.rich, Ktrop$gamma.stat, pch=17, col=Kcol, cex = cexpt)
  points(Ktrop.slice$clade.rich, Ktrop.slice$gamma.stat, pch = 15, col = Kcol.slice, cex = cexpt)
  
  #splines
  points(smooth.spline(Ktrop.slice$clade.rich,Ktrop.slice$gamma.stat,df=5),type='l',col= Kline.slice,lwd=4)
  points(smooth.spline(Ttrop$clade.rich,Ttrop$gamma.stat,df=5),type='l',col= Tline,lwd=4)
  points(smooth.spline(Ktrop$clade.rich,Ktrop$gamma.stat,df=5),type='l',col= Kline,lwd=4)
  
  legend('bottomleft',c('Time','Pre-Equilibrium','Energy Gradient'),pch=c(16,15,17), col=c(Tcol,Kcol.slice,Kcol))
  
  plot(Ttemp$clade.rich, Ttemp$gamma.stat, ylim=range(c(Ktrop.slice$gamma.stat, Ttrop$gamma.stat, Ktrop$gamma.stat)),
       xlim = c(0,800),pch=16, col= Tcol, xlab="Clade Richness",ylab="Gamma", las=1, 
       main='Temperate Origin', cex = cexpt)
  
  points(Ktemp$clade.rich, Ktemp$gamma.stat, pch=17, col= Kcol, cex = cexpt)
  points(Ktemp.slice$clade.rich, Ktemp.slice$gamma.stat, pch = 15, col = Kcol.slice, cex = cexpt)
  
  #splines
  points(smooth.spline(Ktemp.slice$clade.rich,Ktemp.slice$gamma.stat,df=5),type='l',col= Kline.slice,lwd=4)
  points(smooth.spline(Ttemp$clade.rich,Ttemp$gamma.stat,df=5),type='l',col= Tline,lwd=4)
  points(smooth.spline(Ktemp$clade.rich,Ktemp$gamma.stat,df=5),type='l',col= Kline,lwd=4)
  
} # end if statement for (not) plotting disturbance

dev.off()
