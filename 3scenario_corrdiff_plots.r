# Plot env-richness and time-richness correlations as function of clade age

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


Dtrop.all = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.trop.csv',sep=''), header=T)
Ttrop.all = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.csv',sep=''), header=T)
Ktrop.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)

Dtemp.all = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.temp.csv',sep=''), header=T)
Ttemp.all = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.temp.csv',sep=''), header=T)
Ktemp.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)

Dtrop = subset(Dtrop.all, clade.richness >= 30 & n.regions >=5)
Ttrop = subset(Ttrop.all, clade.richness >= 30 & n.regions >=5)
Ktrop = subset(Ktrop.all, clade.richness >= 30 & n.regions >=5)
Dtemp = subset(Dtemp.all, clade.richness >= 30 & n.regions >=5)
Ttemp = subset(Ttemp.all, clade.richness >= 30 & n.regions >=5)
Ktemp = subset(Ktemp.all, clade.richness >= 30 & n.regions >=5)

#Calculate difference between time-richness correlation and environment-richness correlation
Dtrop$time.env.r.diff = Dtrop$r.time.rich - Dtrop$r.env.rich
Ttrop$time.env.r.diff = Ttrop$r.time.rich - Ttrop$r.env.rich
Ktrop$time.env.r.diff = Ktrop$r.time.rich - Ktrop$r.env.rich
Dtemp$time.env.r.diff = Dtemp$r.time.rich - Dtemp$r.env.rich
Ttemp$time.env.r.diff = Ttemp$r.time.rich - Ttemp$r.env.rich
Ktemp$time.env.r.diff = Ktemp$r.time.rich - Ktemp$r.env.rich


pdf(paste(analysis_dir,'/summaryplots/corrDiff__vs_cladeAgeRich',Sys.Date(),'.pdf',sep=''), height=6, width = 8)
par(mfrow=c(2,3), mar = c(5,2,2,1), oma=c(1,5,1,1), mgp = c(4,1,0))
cexpts = 0.5
cexaxis = 1.5
cexlab = 2
cexmain = 2
#Vs clade origin time
plot(log10(Ttrop$clade.origin.time), Ttrop$time.env.r.diff, pch=16, col='red',ylim=c(-2,2),xlab="",
     ylab= "",cex=cexpts, main='Time', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain)
points(log10(Ttemp$clade.origin.time), Ttemp$time.env.r.diff, col='blue',cex=cexpts)
abline(h=0,lty='dashed')
legend('bottomleft',c('temperate origin','tropical origin'), pch = c(1,16), col = c('blue','red'), cex=1.5)

plot(log10(Dtrop$clade.origin.time), Dtrop$time.env.r.diff, pch=16, col='red',ylim=c(-2,2),
     xlab="", ylab= "",cex=cexpts, main='Disturbance', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain)
points(log10(Dtemp$clade.origin.time), Dtemp$time.env.r.diff, col='blue',cex=cexpts)
abline(h=0,lty='dashed')
mtext(expression(paste(plain(log)[10]," Time of Clade Origin")), 1, line=3.5, cex = 1.5)

plot(log10(Ktrop$clade.origin.time), Ktrop$time.env.r.diff, pch=16, col='red',ylim=c(-2,2),xlab="",
     ylab= "",cex=cexpts, main='Energy Gradient', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain)
points(log10(Ktemp$clade.origin.time), Ktemp$time.env.r.diff, col='blue',cex=cexpts)
abline(h=0,lty='dashed')

# Vs clade richness
plot(log10(Ttrop$clade.richness), Ttrop$time.env.r.diff, pch=16, col='red',ylim=c(-2,2),xlab="",
     ylab= "",cex=cexpts, main='', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain)
points(log10(Ttemp$clade.richness), Ttemp$time.env.r.diff, col='blue',cex=cexpts)
abline(h=0,lty='dashed')

plot(log10(Dtrop$clade.richness), Dtrop$time.env.r.diff, pch=16, col='red',ylim=c(-2,2),
     xlab="", ylab= "",cex=cexpts, main='', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain)
points(log10(Dtemp$clade.richness), Dtemp$time.env.r.diff, col='blue',cex=cexpts)
abline(h=0,lty='dashed')
mtext(expression(paste(plain(log)[10]," Clade Richness")), 1, line=3.5, cex = 1.5)

plot(log10(Ktrop$clade.richness), Ktrop$time.env.r.diff, pch=16, col='red',ylim=c(-2,2),xlab="",
     ylab= "",cex=cexpts, main='', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain)
points(log10(Ktemp$clade.richness), Ktemp$time.env.r.diff, col='blue',cex=cexpts)
abline(h=0,lty='dashed')

mtext(expression(paste(italic(r)[time-richness]," - ",italic(r)[environment-richness])), 2, outer=T, line=2, cex = cexlab)

dev.off()
