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

Dcol = 'gold2'
Tcol = 'olivedrab3'
Kcol = 'mediumorchid2'
Dline = 'goldenrod3'
Tline = 'darkgreen'
Kline = 'mediumorchid4'

pdf(paste(analysis_dir,'/summaryplots/3scenarios_gamma_logcladerich_',Sys.Date(),'.pdf',sep=''),height=5, width=10)
par(mfrow=c(1,2), mar = c(4,4,2,1))
cexpt = .8
plot(log10(Ttrop$clade.rich), Ttrop$gamma.stat, ylim=range(c(Dtrop$gamma.stat, Ttrop$gamma.stat, Ktemp$gamma.stat)),
     xlim = c(.9,4),pch=16, col=Tcol,xlab="Clade Richness",ylab="Gamma", las=1, xaxt="n",
     main='Tropical Origin', cex = cexpt)
axis(1,1:4)
points(log10(Dtrop$clade.rich), Dtrop$gamma.stat, pch = 15, col = Dcol, cex = cexpt)
points(log10(Ktrop$clade.rich), Ktrop$gamma.stat, pch=17, col=Kcol, cex = cexpt)
abline(h=0, lty='dashed')

#splines
points(smooth.spline(log10(Dtrop$clade.rich),Dtrop$gamma.stat,df=5),type='l',col= Dline,lwd=4)
points(smooth.spline(log10(Ttrop$clade.rich),Ttrop$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(log10(Ktrop$clade.rich),Ktrop$gamma.stat,df=5),type='l',col= Kline,lwd=4)

legend('bottomleft',c('Time','Disturbance','Energy Gradient'),pch=c(16,15,17), col=c(Tcol,Dcol,Kcol))

plot(log10(Ttemp$clade.rich), Ttemp$gamma.stat, ylim=range(c(Dtrop$gamma.stat, Ttrop$gamma.stat, Ktemp$gamma.stat)),
     xlim = c(0.9,4),pch=16, col= Tcol, xlab="Clade Richness",ylab="Gamma", las=1, xaxt="n",
     main='Temperate Origin', cex = cexpt)
axis(1,1:4)
points(log10(Dtemp$clade.rich), Dtemp$gamma.stat, pch = 15, col = Dcol, cex = cexpt)
points(log10(Ktemp$clade.rich), Ktemp$gamma.stat, pch=17, col= Kcol, cex = cexpt)
abline(h=0, lty='dashed')

#splines
points(smooth.spline(log10(Dtemp$clade.rich),Dtemp$gamma.stat,df=5),type='l',col= Dline,lwd=4)
points(smooth.spline(log10(Ttemp$clade.rich),Ttemp$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(log10(Ktemp$clade.rich),Ktemp$gamma.stat,df=5),type='l',col= Kline,lwd=4)

dev.off()
