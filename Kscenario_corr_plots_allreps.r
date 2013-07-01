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


Ktrop.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)

Ktemp.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)

Ktrop = subset(Ktrop.all, clade.richness >= 30 & n.regions >=5)
Ktemp = subset(Ktemp.all, clade.richness >= 30 & n.regions >=5)

pdf(paste(analysis_dir,'/summaryplots/Kcorr_vs_cladeAgeRich',Sys.Date(),'.pdf',sep=''), height=6, width = 12)
par(mfrow=c(1,2), mar = c(5,5.5,1,1), oma=c(1,1,1,1), mgp = c(4,1,0))
cexpts = 0.5
cexaxis = 1.5
cexlab = 1.5
cexmain = 2
#Vs clade origin time
plot(log10(Ktrop$clade.origin.time), Ktrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab= "Environment-richness correlation",cex=cexpts, main='', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain, las=1)
points(log10(Ktemp$clade.origin.time), Ktemp$r.env.rich, col='blue',cex=cexpts)
abline(h=0,lty='dashed')
legend('bottomleft',c('temperate origin','tropical origin'), pch = c(1,16), col = c('blue','red'), cex=1.5)
mtext(expression(paste(plain(log)[10]," Time of Clade Origin")),1,outer=F, line=3, cex = cexlab)
mtext("(a)",3,at=-1,cex=2)

#Vs clade richness
plot(log10(Ktrop$clade.richness), Ktrop$r.env.rich, pch=16, col='red',ylim=c(-1,1),xlab="",
     ylab= "",cex=cexpts, main='', cex.lab = cexlab, cex.axis = cexaxis, cex.main=cexmain, las=1)
points(log10(Ktemp$clade.richness), Ktemp$r.env.rich, col='blue',cex=cexpts)
abline(h=0,lty='dashed')
mtext(expression(paste(plain(log)[10]," Clade Richness")),1,outer=F, line=3, cex = cexlab)
mtext("(b)",3,at=1,cex=2)
dev.off()


