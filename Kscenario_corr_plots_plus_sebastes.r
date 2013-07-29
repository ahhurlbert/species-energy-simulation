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
  sim_dir = "C:/SENCoutput/senc_reps_analysis"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
  repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"
}

if (Allen == 0) {
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204"
  repo_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation"  
}

#Simulation data
Ktrop.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.trop.csv',sep=''), header=T)
Ktemp.all = read.csv(paste(sim_dir,'/SENC_Stats_K.sims.temp.csv',sep=''), header=T)

Ktrop = subset(Ktrop.all, clade.richness >= 30 & n.regions >=5)
Ktemp = subset(Ktemp.all, clade.richness >= 30 & n.regions >=5)

#Sebastes data
sebastes = read.csv(paste(repo_dir,'/sebastes_data_for_allen.csv',sep=''),header=T)
phy = read.tree(paste(repo_dir,'/Sebastes_tree_Ingram2011PRSB.phy',sep=''))
#Drop non-NEP species (with no latitude data)
nonNEPsp = as.character(sebastes[is.na(sebastes$min_latitude), 'X'])
NEPphy = drop.tip(phy,nonNEPsp)
lat = min(sebastes$min_latitude, na.rm=T):max(sebastes$max_latitude, na.rm=T)

min.num.spp = 5

lat.corr.output = c()
for (c in (NEPphy$Nnode+2):max(NEPphy$edge)) {
  
  #pull out list of species names belonging to each subclade
  sub.clade = clade.members(c, NEPphy, tip.labels=T)
  sub.populations = subset(sebastes, X %in% sub.clade);
  
  sub.richness = sapply(lat, function(x) nrow(sub.populations[sub.populations$min_latitude <= x & sub.populations$max_latitude >= x, ]))
  if(length(sub.clade) >= min.num.spp) {
    lat.corr = cor(lat[sub.richness>0], sub.richness[sub.richness>0])
    lat.corr2 = cor(lat[sub.richness>0 & lat >= 34], sub.richness[sub.richness > 0 & lat >=34])
    lat.corr3 = cor(lat[sub.richness>0 & lat < 34], sub.richness[sub.richness > 0 & lat < 34])
    lat.corr.output = rbind(lat.corr.output, c(c, length(sub.clade), lat.corr, lat.corr2, lat.corr3))
  }
}
lat.corr.output = data.frame(lat.corr.output)
names(lat.corr.output) = c('cladeID','clade.richness','r.lat.rich','r.lat.rich.gte34','r.lat.rich.lt34')  

# Example clade richness vs opportunity patterns
source('opportunity_gradient.r')
sim.for.plot = 3325
stats = read.csv(paste('C:/SENCoutput/SENC_stats_sim',sim.for.plot,'.csv',sep=''), header=T)

c3157 = reg.opp.grad(cladeid = 3157, sim = sim.for.plot, plot=F, output=T)
c2486 = reg.opp.grad(cladeid = 2486, sim = sim.for.plot, plot=F, output=T)
col.3157 = 'gray80' #colors()[577]
col.2486 = 'black' #colors()[149]

#Plot
pdf(paste(analysis_dir,'/summaryplots/latcorr_subclades_vs_cladeAgeRich_',Sys.Date(),'.pdf',sep=''), height=10, width = 14)
par(mfcol=c(2,2), mar = c(7, 8, 3, 4), oma=c(1, 1, 1, 4), mgp = c(5.5, 1.5, 0))
cexpts = 2
cexpts.seb = 3
cexaxis = 2
cexlab = 2.75
cexleg = 2
cexabc = 2.5
pch.temp = 18

#Vs clade origin time (negative the r.env.rich is equal to the richness-latitude correlation)
# plotting complement of x-axis (hence, difference from max value)
plot(max(log10(Ktrop$clade.origin.time)) - log10(Ktrop$clade.origin.time), -Ktrop$r.env.rich, xaxt = "n",
     xlab = expression(paste(plain(log)[10]," Clade origin time")), pch=16, col='red',ylim=c(-1,1),
     ylab = expression(italic(r)[latitude-richness]), cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las=1)
points(max(log10(Ktrop$clade.origin.time)) - log10(Ktemp$clade.origin.time), -Ktemp$r.env.rich, 
           col = 'blue', pch = pch.temp, cex = cexpts)
points(max(log10(Ktrop$clade.origin.time)) - log10(stats$clade.origin.time[stats$clade.id == 3157]),
       -stats$r.env.rich[stats$clade.id == 3157], pch = 17, col = col.3157, cex = 1.5*cexpts)
points(max(log10(Ktrop$clade.origin.time)) - log10(stats$clade.origin.time[stats$clade.id == 2486]),
       -stats$r.env.rich[stats$clade.id == 2486], pch = 17, col = col.2486, cex = 1.5*cexpts)
axis(1, at = 0:4, labels = 4:0, cex.axis = cexaxis)
mtext(c("recent","old"), 1, at = c(0,4), line = 3, cex = 1.5)
abline(h = 0, lty = 'dashed')
legend('topright',c('temperate origin','tropical origin'), pch = c(pch.temp, 16), col = c('blue','red'), cex = cexleg)
mtext("(a)", 2, at = 1.3, cex = cexabc, outer = F, las = 1, line = 5)

#Vs clade richness
plot(log10(Ktrop$clade.richness), -Ktrop$r.env.rich, pch = 16, col = 'red', ylim = c(-1,1),
     xlab = expression(paste(plain(log)[10]," Clade richness")), ylab = expression(italic(r)[latitude-richness]),
     cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las = 1)
points(log10(Ktemp$clade.richness), -Ktemp$r.env.rich, col = 'blue', pch = pch.temp, cex = cexpts)
points(log10(stats$clade.richness[stats$clade.id == 3157]),
       -stats$r.env.rich[stats$clade.id == 3157], pch = 17, col = col.3157, cex = 2*cexpts)
points(log10(stats$clade.richness[stats$clade.id == 2486]),
       -stats$r.env.rich[stats$clade.id == 2486], pch = 17, col = col.2486, cex = 2*cexpts)
abline(h = 0,lty = 'dashed')
mtext("(b)", 2, outer=F, at = 1.3, cex = cexabc, las = 1, line = 5)


#Plot example richness and opportunity gradients for two examples
lat.bin = 11 - c3157$region
eq.rich.grad = lm(c3157$spp.rich ~ lat.bin)

plot(11 - c3157$region,c3157$spp.rich, type='n', xlim=c(1,10), ylim = c(0,500),
     xlab = 'Latitude', ylab = 'Species richness', xaxt = "n", cex.lab = cexlab, cex.axis = cexaxis, las=1)
lines(10:1, predict(eq.rich.grad, newdata = data.frame(1:10)), lwd = 2)
points(11 - c2486$region, c2486$clade.rich, type = 'l', col = col.2486, lwd = 5)
arrows(11 - c2486$reg.of.origin, -0.25*max(c2486$spp.rich), 11 - c2486$reg.of.origin, -20, col = col.2486, lwd = 3, xpd = T)
points(11 - c3157$region, c3157$clade.rich, type = 'l', col = col.3157, lwd = 5)
arrows(11 - c3157$reg.of.origin, -0.25*max(c3157$spp.rich), 11 - c3157$reg.of.origin, -20, col = col.3157, lwd = 3, xpd = T)
mtext("(c)", 2, at = 570, outer = F, cex = cexabc, las = 1, line = 5)
par(new=T)
ymax = 1.2
plot(11-c3157$region, c3157$opp.frac, type = 'l', xlim = c(1,10), col=col.3157, yaxt="n", xaxt = "n",
     xlab="", ylab="", ylim=c(0, ymax), lwd = 5, lty = 'dashed')
points(11-c2486$region, c2486$opp.frac, type = 'l', xlim = c(1,10), col=col.2486, 
     ylim=c(0, ymax), lwd = 5, lty = 'dashed')
axis(4, at = seq(0, 1, by = .2), cex.axis = cexaxis, las = 1)
axis(1, at = 1:10, labels = F)
mtext("Opportunity", 4, outer=T, cex = 2.25, adj = 0.84, line = 2)
mtext("Tropics",1,adj=.05,line=1.5, cex=1.5)
mtext("Temperate",1,adj=.95,line=1.5, cex=1.5)


#Sebastes vs clade richness
plot(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich, ylim = c(-1,1), las=1, cex.axis = cexaxis, cex.lab = cexlab,
     xlab = expression(paste(plain(log)[10]," Clade richness")), ylab = expression(italic(r)[latitude-richness]), pch = 16, cex = cexpts.seb)
points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.gte34, pch = 1, cex = cexpts.seb)
#points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.lt34, pch=2, cex = cexpts)
abline(h=0,lty='dashed')
legend("topright", c('Entire gradient','North of 34N'), pch = c(16, 1), cex = cexleg)
mtext("(d)", 2, outer = F, at = 1.3, cex = cexabc, las = 1, line = 5)
dev.off()








###############
# Clade richness vs time of clade origin, to help interpret above patterns

pdf(paste(analysis_dir,'/summaryplots/cladeRich_vs_origintime_',Sys.Date(),'.pdf',sep=''), height=6, width = 7)
par(mar=c(5,5,1,1))
plot(log10(Ktrop$clade.origin.time),log10(Ktrop$clade.richness), pch=16, col='red', cex.lab = 2,
     xlab = expression(paste(plain(log)[10]," Clade origin time")), 
     ylab = expression(paste(plain(log)[10]," Clade richness")))
points(log10(Ktemp$clade.origin.time),log10(Ktemp$clade.richness),pch=18,col='blue')
dev.off()
