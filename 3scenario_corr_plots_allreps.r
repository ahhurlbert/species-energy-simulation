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


Dtrop = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.trop.csv',sep=''), header=T)
Ttrop = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.csv',sep=''), header=T)
Ktrop = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)

Dtemp = read.csv(paste(sim_dir,'/SENC_Stats_D.sims.temp.csv',sep=''), header=T)
Ttemp = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.temp.csv',sep=''), header=T)
Ktemp = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)


pdf(paste(analysis_dir,'/summaryplots/corr_cladeorigin_plots_',Sys.Date(),'.pdf',sep=''), height=6, width = 8)
par(mfrow=c(2,3), mar = c(2,2,2,1), oma=c(4,5,1,1))
cexpts = 0.5
cexlab = 1.5
cexmain = 2
#environment-richness plots
plot(log10(Ttrop$clade.origin.time), Ttrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab= "",cex=cexpts, main='Time', cex.lab=2, cex.main=cexmain)
points(log10(Ttemp$clade.origin.time), Ttemp$r.env.rich, col='blue',cex=cexpts)
mtext(expression(italic(r)[environment-richness]), 2, at = 0.75, outer=T, line=2, cex = cexlab)

plot(log10(Dtrop$clade.origin.time), Dtrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts, main='Disturbance', cex.main=cexmain)
points(log10(Dtemp$clade.origin.time), Dtemp$r.env.rich, col='blue',cex=cexpts)

plot(log10(Ktrop$clade.origin.time), Ktrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts, main='Energy Gradient', cex.main=cexmain)
points(log10(Ktemp$clade.origin.time), Ktemp$r.env.rich, col='blue',cex=cexpts)

# time-richness plots
plot(log10(Ttrop$clade.origin.time), Ttrop$r.time.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab= "",cex=cexpts, cex.lab=2)
points(log10(Ttemp$clade.origin.time), Ttemp$r.time.rich, col='blue',cex=cexpts)
mtext(expression(italic(r)[time-richness]), 2, at = 0.25, outer=T, line=2, cex = cexlab)
legend('bottomleft',c('temperate origin','tropical origin'), pch = c(1,16), col = c('blue','red'), cex=1)

plot(log10(Dtrop$clade.origin.time), Dtrop$r.time.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts)
points(log10(Dtemp$clade.origin.time), Dtemp$r.time.rich, col='blue',cex=cexpts)

plot(log10(Ktrop$clade.origin.time), Ktrop$r.time.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts)
points(log10(Ktemp$clade.origin.time), Ktemp$r.time.rich, col='blue',cex=cexpts)

mtext(expression(paste(plain(log)[10]," Time of Clade Origin")),1,outer=T, line=2, cex = cexlab)
dev.off()

########################

# Vs clade richness instead of clade age

pdf(paste(analysis_dir,'/summaryplots/corr_claderichness_plots_',Sys.Date(),'.pdf',sep=''), height=6, width = 8)
par(mfrow=c(2,3), mar = c(2,2,2,1), oma=c(4,5,1,1))
cexpts = 0.5
cexlab = 1.5
cexmain = 2
#environment-richness plots
plot(log10(Ttrop$clade.richness), Ttrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab= "",cex=cexpts, main='Time', cex.lab=2, cex.main=cexmain)
points(log10(Ttemp$clade.richness), Ttemp$r.env.rich, col='blue',cex=cexpts)
mtext(expression(italic(r)[environment-richness]), 2, at = 0.75, outer=T, line=2, cex = cexlab)

plot(log10(Dtrop$clade.richness), Dtrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts, main='Disturbance', cex.main=cexmain)
points(log10(Dtemp$clade.richness), Dtemp$r.env.rich, col='blue',cex=cexpts)

plot(log10(Ktrop$clade.richness), Ktrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts, main='Energy Gradient', cex.main=cexmain)
points(log10(Ktemp$clade.richness), Ktemp$r.env.rich, col='blue',cex=cexpts)

# time-richness plots
plot(log10(Ttrop$clade.richness), Ttrop$r.time.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab= "",cex=cexpts, cex.lab=2)
points(log10(Ttemp$clade.richness), Ttemp$r.time.rich, col='blue',cex=cexpts)
mtext(expression(italic(r)[time-richness]), 2, at = 0.25, outer=T, line=2, cex = cexlab)
legend('bottomleft',c('temperate origin','tropical origin'), pch = c(1,16), col = c('blue','red'), cex=1)

plot(log10(Dtrop$clade.richness), Dtrop$r.time.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts)
points(log10(Dtemp$clade.richness), Dtemp$r.time.rich, col='blue',cex=cexpts)

plot(log10(Ktrop$clade.richness), Ktrop$r.time.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab="",cex=cexpts)
points(log10(Ktemp$clade.richness), Ktemp$r.time.rich, col='blue',cex=cexpts)

mtext(expression(paste(plain(log)[10]," Time of Clade Origin")),1,outer=T, line=2, cex = cexlab)
dev.off()

