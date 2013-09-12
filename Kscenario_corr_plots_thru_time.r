# Figure 4: Plot latitude-richness correlation versus (a) time of clade origin and
# (b) clade richness for all simulation subclades with at least 30 species and
# spanning at least 5 regions. (c) Plot same thing for empirical Sebastes data
# for the entire northeastern Pacific gradient (23-66N) as well as the gradient
# north of Point Conception (34-66N).

# libraries
require(ape)
require(caper)


Allen = 1;

if (Allen ==1) {
  sim_dir = "C:/SENCoutput/longtime_reps"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
  repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"
}

if (Allen == 0) {
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  repo_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation"  
}

#Simulation data
trop.sims = c(4065:4072,4074)
temp.sims = 4075:4084

Ktrop.all = c()
for (i in trop.sims) {
  tmp = read.csv(paste(sim_dir,'/SENC_Stats_sim', i, '_specific_times.csv',sep=''), header=T)
  Ktrop.all = rbind(Ktrop.all, tmp)
}  
Ktemp.all = c()
for (i in temp.sims) {
  tmp = read.csv(paste(sim_dir,'/SENC_Stats_sim', i, '_specific_times.csv',sep=''), header=T)
  Ktemp.all = rbind(Ktemp.all, tmp)
}  

Ktrop = subset(Ktrop.all, clade.richness >= 30 & n.regions >=5)
Ktemp = subset(Ktemp.all, clade.richness >= 30 & n.regions >=5)

times = c(5000, 10000, 15000, 30000, 60000, 100000)

#Plot
pdf(paste(analysis_dir,'/summaryplots/latcorr_subclades_vs_cladeAgeRich_thruTime_',Sys.Date(),'_log.pdf',sep=''), height=10, width = 25)
par(mfcol=c(2,length(times)), mar = c(7, 8, 3, 4), oma=c(1, 1, 1, 4), mgp = c(5.5, 1.5, 0))
cexpts = 2
cexpts.seb = 3
cexaxis = 2
cexlab = 2.75
cexleg = 2
cexabc = 2.5
pch.temp = 18
linewd = 7
arrowwd = 4

for (t in times) {
  Ktrop.sub = subset(Ktrop, time == t)
  Ktemp.sub = subset(Ktemp, time == t)
  #Vs clade origin time (negative the r.env.rich is equal to the richness-latitude correlation)
  # plotting complement of x-axis (hence, difference from max value)
  plot(max(log10(Ktrop.sub$clade.origin.time)) - log10(Ktrop.sub$clade.origin.time), -Ktrop.sub$r.env.rich, xaxt = "n",
       xlab = expression(paste(plain(log)[10]," Clade origin time")), pch=16, col='red',ylim=c(-1,1),
       ylab = "", cex = cexpts, main = paste(t/1000,'k timesteps'), cex.lab = cexlab, cex.axis = cexaxis, las=1, cex.main = 2)
  points(max(log10(Ktrop.sub$clade.origin.time)) - log10(Ktemp.sub$clade.origin.time), -Ktemp.sub$r.env.rich, 
             col = 'blue', pch = pch.temp, cex = cexpts)
  axis(1, at = 0:4, labels = 4:0, cex.axis = cexaxis)
  mtext(c("recent","old"), 1, at = c(0,4), line = 3, cex = 1.5)
  mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)
  abline(h = 0, lty = 'dashed')
  legend('topright',c('temperate origin','tropical origin'), pch = c(pch.temp, 16), col = c('blue','red'), cex = cexleg)

  #Vs clade richness
  plot(log10(Ktrop.sub$clade.richness), -Ktrop.sub$r.env.rich, pch = 16, col = 'red', ylim = c(-1.15,1),
       xlab = expression(paste(plain(log)[10]," Clade richness")), ylab = "",
       cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las = 1)
  points(log10(Ktemp.sub$clade.richness), -Ktemp.sub$r.env.rich, col = 'blue', pch = pch.temp, cex = cexpts)
  abline(h = 0,lty = 'dashed')
  mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)

  #extra tick marks showing % of max richness
  pcts = c(.75,.25,.1)
  axis(1,at=log10(pcts*max(Ktrop.sub$clade.richness)), labels=F,tck= .01)
  text(log10(pcts*max(Ktrop.sub$clade.richness)), rep(-1.12,length(pcts)), paste(pcts*100,"%",sep=""), cex = 1.5) 
}  
dev.off()

