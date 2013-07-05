# Plot phylogenetic gamma statistic as a function of clade richness for the 3 different scenarios
# Energy Gradient, no zero-sum, and Energy Gradient pre-equilibrium
# a) tropical origin, b) temperate origin, c) comparison of simulation data to gamma for
# Sebastes rockfish phylogeny

require(ape)

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

# Empirical data for Sebastes
phy = read.tree(paste(repo_dir,'/Sebastes_tree_Ingram2011PRSB.phy',sep=''))
sebastes = read.csv(paste(repo_dir,'/sebastes_data_for_allen.csv',sep=''),header=T)
#Drop non-NEP species (with no latitude data)
nonNEPsp = as.character(sebastes[is.na(sebastes$min_latitude), 'X'])
NEPphy = drop.tip(phy,nonNEPsp)

# Simulation data
Ttrop = read.csv(paste(sim_dir,'/SENC_Stats_T.sims.trop.csv',sep=''), header=T)
Ktrop = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)
Ktrop.slice = read.csv(paste(sim_dir,'/SENC_Stats_K.slice.sims.trop.csv',sep=''), header=T)

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

pdf(paste(analysis_dir,'/summaryplots/3scenarios_gamma_sim_and_empirical_',Sys.Date(),'.pdf',sep=''), height = 15, width = 7)
par(mfrow=c(3,1), mar = c(6.5,8,2,1), mgp = c(5.5, 1.5, 0), oma = c(1, 5, 1, 1))

#cexaxis = 2

cexpt = 2
cexlab = 3
cexaxis = 2.25
cexlegend = 2.5
cexabc = 2.5

plot(log10(Ttrop$clade.rich), Ttrop$gamma.stat, xlim = c(0.9,4), ylim = c(-20,5), pch=16, col=Tcol,
     xlab=expression(paste(plain(log)[10]," Clade Richness")), ylab = "Gamma", las = 1, 
     xaxt="n", main='', cex = cexpt, cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
axis(1,1:4, cex.axis=cexaxis)
points(log10(Ktrop$clade.rich), Ktrop$gamma.stat, pch=17, col=Kcol, cex = cexpt)
points(log10(Ktrop.slice$clade.rich), Ktrop.slice$gamma.stat, pch = 15, col = Kcol.slice, cex = cexpt)
  
#splines
points(smooth.spline(log10(Ktrop.slice$clade.rich),Ktrop.slice$gamma.stat,df=5),type='l',col= Kline.slice,lwd=4)
points(smooth.spline(log10(Ttrop$clade.rich),Ttrop$gamma.stat,df=5),type='l',col= Tline,lwd=4)
points(smooth.spline(log10(Ktrop$clade.rich),Ktrop$gamma.stat,df=5),type='l',col= Kline,lwd=4)
  
legend('bottomleft',c('Time','Pre-Equilibrium','Energy Gradient'),
       pch = c(16,15,17), col = c(Tcol, Kcol.slice, Kcol), bty = "n", cex = cexlegend)
mtext("(a)", 2, at = 1, cex = cexabc, outer = T, las = 1, line = .5)

#Insets with cartoon phylogenies
posgamma = read.tree(paste(repo_dir,'/posgamma.tre',sep=''))
neggamma = read.tree(paste(repo_dir,'/neggamma.tre',sep=''))
par(new = T, fig = c(.65,1,.85,.99))
plot(posgamma, edge.width=3)
par(new = T, fig = c(.65,1,.68,.82))
plot(neggamma, edge.width=3)
  
par(new=T)
par(fig = c(0, 1, 0.33, .67))
  
plot(log10(Ttemp$clade.rich), Ttemp$gamma.stat, ylim = c(-20,5), xlim = c(0.9,4), pch=16, col= Tcol, 
     xlab=expression(paste(plain(log)[10]," Clade Richness")), ylab = "Gamma", las = 1, xaxt = "n",
     main='', cex = cexpt, cex.lab = cexlab, cex.axis = cexaxis, cex.main = cexmain)
axis(1, 1:4, cex.axis = cexaxis)
points(log10(Ktemp$clade.rich), Ktemp$gamma.stat, pch = 17, col = Kcol, cex = cexpt)
points(log10(Ktemp.slice$clade.rich), Ktemp.slice$gamma.stat, pch = 15, col = Kcol.slice, cex = cexpt)
  
#splines
points(smooth.spline(log10(Ktemp.slice$clade.rich), Ktemp.slice$gamma.stat, df = 5), type = 'l', col = Kline.slice, lwd=4)
points(smooth.spline(log10(Ttemp$clade.rich), Ttemp$gamma.stat, df = 5), type = 'l', col = Tline, lwd=4)
points(smooth.spline(log10(Ktemp$clade.rich), Ktemp$gamma.stat, df = 5), type = 'l', col = Kline, lwd=4)
  
mtext("(b)", 2, at = 0.67, cex = cexabc, outer = T, las = 1, line = .5)

#Empirical data
#Sebastes phylogeny has 99 species (only 66 in NEP), so pull out clades for each scenario of roughly the same size
rich=66
  
Ttrop99 = subset(Ttrop, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktrop99 = subset(Ktrop, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktrop.slice99 = subset(Ktrop.slice, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktemp99 = subset(Ktemp, clade.richness > .9*rich & clade.richness < 1.1*rich)
Ktemp.slice99 = subset(Ktemp.slice, clade.richness > .9*rich & clade.richness < 1.1*rich)
  
par(new=T)
par(fig = c(0, 1, 0, 0.33))

plot(density(Ttrop99$gamma.stat), col=Tcol, main="", xlab="Gamma", lwd=5, las=1,
     xlim = c(-8,2), cex.lab = cexlab, cex.axis = cexaxis)
points(density(Ktrop99$gamma.stat), type='l',col=Kcol, lty='dashed',lwd=5)
points(density(Ktemp99$gamma.stat), type='l',col=Kcol, lwd=5)
points(density(Ktrop.slice99$gamma.stat), type= 'l', col=Kcol.slice, lty='dashed', lwd=5)
points(density(Ktemp.slice99$gamma.stat), type = 'l', col=Kcol.slice, lwd=5)
abline(v = gammaStat(NEPphy), lwd=3)
legend("topleft",c('tropical', 'temperate'),
       col = 'gray50', lty = c('dashed', 'solid'), lwd=4, bty = "n", cex = cexlegend)

mtext("(c)", 2, at = 0.33, cex = cexabc, outer = T, las = 1, line = .5)
dev.off()
