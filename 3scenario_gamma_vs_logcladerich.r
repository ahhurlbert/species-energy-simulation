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


Dtrop.all = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.trop.all.csv',sep=''), header=T)
Ttrop.all = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.all.csv',sep=''), header=T)
Ktrop.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.all.csv',sep=''), header=T)

Dtemp.all = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.temp.all.csv',sep=''), header=T)
Ttemp.all = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.temp.all.csv',sep=''), header=T)
Ktemp.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.all.csv',sep=''), header=T)

Dcol = 'gold2'
Tcol = 'olivedrab3'
Kcol = 'mediumorchid2'
Dline = 'goldenrod3'
Tline = 'darkgreen'
Kline = 'mediumorchid4'

pdf(paste(analysis_dir,'/summaryplots/3scenarios_gamma_logcladerich_',Sys.Date(),'.pdf',sep=''),height=5, width=10)
par(mfrow=c(1,2), mar = c(4,5,2,1))
cexpt = .8
cexlab = 1.5
cexaxis = 1.25
cexmain = 1.25
plot(log10(Ttrop.all$clade.rich), Ttrop.all$gamma.stat, ylim=range(c(Dtrop.all$gamma.stat, Ttrop.all$gamma.stat, Ktemp.all$gamma.stat)),
     xlim = c(.9,4),pch=16, col=Tcol,xlab=expression(paste(plain(log)[10]," Clade Richness")),ylab="Gamma", las=1, xaxt="n",
     main='Tropical Origin', cex = cexpt, cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
axis(1,1:4, cex.axis = cexaxis)
points(log10(Dtrop.all$clade.rich), Dtrop.all$gamma.stat, pch = 15, col = Dcol, cex = cexpt)
points(log10(Ktrop.all$clade.rich), Ktrop.all$gamma.stat, pch=17, col=Kcol, cex = cexpt)
abline(h=0, lty='dashed')

#splines
points(smooth.spline(log10(Dtrop.all$clade.rich),Dtrop.all$gamma.stat,df=5),type='l',col= Dline,lwd=4)
points(smooth.spline(log10(Ttrop.all$clade.rich),Ttrop.all$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(log10(Ktrop.all$clade.rich),Ktrop.all$gamma.stat,df=5),type='l',col= Kline,lwd=4)

legend('bottomleft',c('Time','Disturbance','Energy Gradient'),pch=c(16,15,17), col=c(Tcol,Dcol,Kcol))

#Insets with cartoon phylogenies
require(ape)
posgamma = read.tree(paste(repo_dir,'/posgamma.tre',sep=''))
neggamma = read.tree(paste(repo_dir,'/neggamma.tre',sep=''))
par(new=T)
par(fig = c(.3,.5,.6,.98))
plot(posgamma, edge.width=2)
par(new=T)
par(fig = c(.3,.5,.02,.4))
plot(neggamma, edge.width=2)

par(new=T)
par(fig = c(.5,1,0,1))

plot(log10(Ttemp.all$clade.rich), Ttemp.all$gamma.stat, ylim=range(c(Dtrop.all$gamma.stat, Ttrop.all$gamma.stat, Ktemp.all$gamma.stat)),
     xlim = c(0.9,4),pch=16, col= Tcol, xlab=expression(paste(plain(log)[10]," Clade Richness")),ylab="Gamma", las=1, xaxt="n",
     main='Temperate Origin', cex = cexpt, cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
axis(1,1:4, cex.axis=cexaxis)
points(log10(Dtemp.all$clade.rich), Dtemp.all$gamma.stat, pch = 15, col = Dcol, cex = cexpt)
points(log10(Ktemp.all$clade.rich), Ktemp.all$gamma.stat, pch=17, col= Kcol, cex = cexpt)
abline(h=0, lty='dashed')

#splines
points(smooth.spline(log10(Dtemp.all$clade.rich),Dtemp.all$gamma.stat,df=5),type='l',col= Dline,lwd=4)
points(smooth.spline(log10(Ttemp.all$clade.rich),Ttemp.all$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(log10(Ktemp.all$clade.rich),Ktemp.all$gamma.stat,df=5),type='l',col= Kline,lwd=4)

dev.off()
