# Code for plotting Figure 3 in Hurlbert & Stegen

# libraries
require(ape)
require(caper)

##########################################
# Specify simulations for plotting
##########################################

#Energy gradient sims
# currently has 100 sims
Ktrop.sims = c(4065:4074, 4185:4274)
Ktemp.sims = c(4075:4084, 4275:4364)

#No energy gradient sims
# currently has 100 sims
Ttrop.sims = 3465:3564
Ttemp.sims = 3565:3664

#########################################


#Simulation data
compile.stats = function(sims) {
  for (sim in sims) {
    compiled = c()
    temp = read.csv(paste('raw_sim_output/Stats_sim', sim, '_all_subclades.csv', sep = ''), header = T)
    compiled = rbind(compiled, temp)
  }
}

Ktrop.all = compile.stats(Ktrop.sims)
Ktemp.all = compile.stats(Ktemp.sims)
Ttrop.all = compile.stats(Ttrop.sims)
Ttemp.all = compile.stats(Ttemp.sims)

Ttrop = subset(Ttrop.all, clade.richness >= 30 & n.regions >=5)
Ktrop = subset(Ktrop.all, clade.richness >= 30 & n.regions >=5)
Ttemp = subset(Ttemp.all, clade.richness >= 30 & n.regions >=5)
Ktemp = subset(Ktemp.all, clade.richness >= 30 & n.regions >=5)

#####################################################################
#Sebastes data
sebastes = read.csv('sebastes_data_for_allen.csv', header = T)
phy = read.tree('Sebastes_tree_Ingram2011PRSB.phy')

#Drop non-NorthEastern Pacific species (with no latitude data)
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
    lat.corr = cor(lat[sub.richness > 0], sub.richness[sub.richness > 0])
    lat.corr2 = cor(lat[sub.richness > 0 & lat >= 34], sub.richness[sub.richness > 0 & lat >= 34])
    lat.corr3 = cor(lat[sub.richness > 0 & lat < 34], sub.richness[sub.richness > 0 & lat < 34])
    lat.corr.output = rbind(lat.corr.output, c(c, length(sub.clade), lat.corr, lat.corr2, lat.corr3))
  }
}
lat.corr.output = data.frame(lat.corr.output)
names(lat.corr.output) = c('cladeID','clade.richness','r.lat.rich','r.lat.rich.gte34','r.lat.rich.lt34')  

####################################################################

# Example clade richness vs opportunity patterns
source('opportunity_gradient.r')
sim.for.plot = 3325
stats = read.csv(paste('analysis_output/SENC_stats_sim',sim.for.plot,'.csv',sep=''), header=T)

c3157 = reg.opp.grad(cladeid = 3157, sim = sim.for.plot, plot=F, output=T)
c2486 = reg.opp.grad(cladeid = 2486, sim = sim.for.plot, plot=F, output=T)
col.3157 = 'gray80' #colors()[577]
col.2486 = 'black' #colors()[149]

#Plot
pdf('analysis_output/latcorr_subclades_vs_cladeAgeRich.pdf', height=15, width = 14)
par(mfrow=c(3,2), mar = c(7, 8, 3, 10), oma=c(1, 1, 1, 1), mgp = c(5.5, 1.5, 0))
cexpts = 2
cexpts.seb = 3
cexaxis = 2
cexlab = 2.75
cexleg = 2
cexabc = 2.5
pch.temp = 18
linewd = 7
arrowwd = 4

#No energy gradient plots
plot(max(log10(Ttrop$clade.origin.time)) - log10(Ttrop$clade.origin.time), -Ttrop$r.env.rich, xaxt = "n",
     xlab = expression(paste(plain(log)[10]," Clade origin time")), pch=16, col='red',ylim=c(-1,1),
     ylab = "", cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las=1)
points(max(log10(Ttrop$clade.origin.time)) - log10(Ttemp$clade.origin.time), -Ttemp$r.env.rich, 
       col = 'blue', pch = pch.temp, cex = cexpts)
axis(1, at = max(log10(Ttemp$clade.origin.time)) - 0:2, labels = 0:2, cex.axis = cexaxis)
mtext(c("recent","old"), 1, at = c(0,2.3), line = 3, cex = 1.5)
mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)
abline(h = 0, lty = 'dashed')
mtext("(a)", 2, at = 1.27, cex = cexabc, outer = F, las = 1, line = 5)
legend(1.2,.7,c('temperate origin','tropical origin'), pch = c(pch.temp, 16), col = c('blue','red'), cex = cexleg)

#Vs clade richness
plot(log10(Ttrop$clade.richness), -Ttrop$r.env.rich, pch = 16, col = 'red', ylim = c(-1.15,1),
     xlab = expression(paste(plain(log)[10]," Clade richness")), ylab = "",
     cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las = 1)
points(log10(Ttemp$clade.richness), -Ttemp$r.env.rich, col = 'blue', pch = pch.temp, cex = cexpts)
abline(h = 0,lty = 'dashed')
mtext("(b)", 2, outer=F, at = 1.3, cex = cexabc, las = 1, line = 5)
mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)

#extra tick marks showing % of max richness
pcts = c(.75, .25, .1, .02)
axis(1,at=log10(pcts*max(Ttrop$clade.richness)), labels=F,tck= .01)
text(log10(pcts*max(Ttrop$clade.richness)), rep(-1.12,length(pcts)), paste(pcts*100,"%",sep=""), cex = 1.75) 

#Energy gradient scenario: Vs clade origin time (negative the r.env.rich is equal to the richness-latitude correlation)
# plotting complement of x-axis (hence, difference from max value)
plot(max(log10(Ktrop$clade.origin.time)) - log10(Ktrop$clade.origin.time), -Ktrop$r.env.rich, xaxt = "n",
     xlab = expression(paste(plain(log)[10]," Clade origin time")), pch=16, col='red',ylim=c(-1,1), xlim = c(.65, 4.75),
     ylab = "", cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las=1)
points(max(log10(Ktrop$clade.origin.time)) - log10(Ktemp$clade.origin.time), -Ktemp$r.env.rich, 
           col = 'blue', pch = pch.temp, cex = cexpts)
points(max(log10(Ktrop$clade.origin.time)) - log10(stats$clade.origin.time[stats$clade.id == 3157]),
       -stats$r.env.rich[stats$clade.id == 3157], pch = 17, col = 'white', cex = 2.5*cexpts)
points(max(log10(Ktrop$clade.origin.time)) - log10(stats$clade.origin.time[stats$clade.id == 2486]),
       -stats$r.env.rich[stats$clade.id == 2486], pch = 17, col = 'white', cex = 2.5*cexpts)
points(max(log10(Ktrop$clade.origin.time)) - log10(stats$clade.origin.time[stats$clade.id == 2486]),
       -stats$r.env.rich[stats$clade.id == 2486], pch = 17, col = col.2486, cex = 2*cexpts)
points(max(log10(Ktrop$clade.origin.time)) - log10(stats$clade.origin.time[stats$clade.id == 3157]),
       -stats$r.env.rich[stats$clade.id == 3157], pch = 17, col = col.3157, cex = 2*cexpts)
axis(1, at = max(log10(Ktrop$clade.origin.time)) - 4:0, labels = 4:0, cex.axis = cexaxis)
mtext(c("recent","old"), 1, at = c(0,4), line = 3, cex = 1.5)
mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)
abline(h = 0, lty = 'dashed')
mtext("(c)", 2, at = 1.27, cex = cexabc, outer = F, las = 1, line = 5)

#Vs clade richness
plot(log10(Ktrop$clade.richness), -Ktrop$r.env.rich, pch = 16, col = 'red', ylim = c(-1.15,1),
     xlab = expression(paste(plain(log)[10]," Clade richness")), ylab = "",
     cex = cexpts, main = '', cex.lab = cexlab, cex.axis = cexaxis, las = 1)
points(log10(Ktemp$clade.richness), -Ktemp$r.env.rich, col = 'blue', pch = pch.temp, cex = cexpts)
points(log10(stats$clade.richness[stats$clade.id == 3157]),
       -stats$r.env.rich[stats$clade.id == 3157], pch = 17, col = 'white', cex = 2.5*cexpts)
points(log10(stats$clade.richness[stats$clade.id == 2486]),
       -stats$r.env.rich[stats$clade.id == 2486], pch = 17, col = 'white', cex = 2.5*cexpts)
points(log10(stats$clade.richness[stats$clade.id == 2486]),
       -stats$r.env.rich[stats$clade.id == 2486], pch = 17, col = col.2486, cex = 2*cexpts)
points(log10(stats$clade.richness[stats$clade.id == 3157]),
       -stats$r.env.rich[stats$clade.id == 3157], pch = 17, col = col.3157, cex = 2*cexpts)
abline(h = 0,lty = 'dashed')
mtext("(d)", 2, outer=F, at = 1.3, cex = cexabc, las = 1, line = 5)
mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)

#extra tick marks showing % of max richness
axis(1,at=log10(pcts*max(Ktrop$clade.richness)), labels=F,tck= .01)
text(log10(pcts*max(Ktrop$clade.richness)), rep(-1.12,length(pcts)), paste(pcts*100,"%",sep=""), cex = 1.75) 

#Plot example richness and opportunity gradients for two examples
lat.bin = 11 - c3157$region
eq.rich.grad = lm(c3157$spp.rich ~ lat.bin)

plot(11 - c3157$region,log10(c3157$spp.rich), type='n', xlim=c(1,10), ylim = c(0,3), yaxt = "n",
     xlab = 'Latitude', ylab = "", xaxt = "n", cex.lab = cexlab, cex.axis = cexaxis, las=1)
axis(2, at = 0:3, cex.axis = cexaxis, las = 1)
axis(2, at = c(0.5, 1.5, 2.5), labels=F)
lines(10:1, log10(predict(eq.rich.grad, newdata = data.frame(1:10))), lwd = 2)
points(11 - c2486$region, log10(c2486$clade.rich), type = 'l', col = col.2486, lwd = linewd)
arrows(11 - c2486$reg.of.origin, -0.6, 11 - c2486$reg.of.origin, -.15, col = col.2486, lwd = arrowwd, xpd = T)
points(11 - c3157$region, log10(c3157$clade.rich), type = 'l', col = col.3157, lwd = linewd)
arrows(11 - c3157$reg.of.origin, -0.6, 11 - c3157$reg.of.origin, -.15, col = col.3157, lwd = arrowwd, xpd = T)
mtext("(e)", 2, at = 3.45, outer = F, cex = cexabc, las = 1, line = 5)
par(new=T)
ymax = 1.2
plot(11-c3157$region, c3157$opp.frac, type = 'l', xlim = c(1,10), col=col.3157, yaxt="n", xaxt = "n",
     xlab="", ylab="", ylim=c(0, ymax), lwd = linewd, lty = 'dashed')
points(11-c2486$region, c2486$opp.frac, type = 'l', xlim = c(1,10), col=col.2486, 
     ylim=c(0, ymax), lwd = linewd, lty = 'dashed')
axis(4, at = seq(0, 1, by = .2), cex.axis = cexaxis, las = 1)
axis(1, at = 1:10, labels = F)
mtext("Opportunity", 4, outer=F, cex = 2.25, at = .5, line = 5.5)
mtext("Tropics",1,adj=.05,line=1.5, cex=1.5)
mtext("Temperate",1,adj=.95,line=1.5, cex=1.5)
mtext(expression(paste(plain(log)[10]," Clade richness")), 2, cex = 2, line = 4.5)

#Sebastes vs clade richness
seb.col1 = 'olivedrab3'
seb.col2 = 'darkgreen'
plot(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich, ylim = c(-1.15,1), las=1, cex.axis = cexaxis, cex.lab = cexlab,
     xlab = expression(paste(plain(log)[10]," Clade richness")), ylab = "", cex = cexpts.seb, col = seb.col1, pch = 16)
points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.gte34, pch = 16, cex = cexpts.seb, col = seb.col2)
#points(log10(lat.corr.output$clade.richness), lat.corr.output$r.lat.rich.lt34, pch=2, cex = cexpts)
abline(h=0,lty='dashed')
legend("topright", c('Entire gradient','North of 34N'), pch = c(16, 16), cex = cexleg, col = c(seb.col1, seb.col2))
mtext("(f)", 2, outer = F, at = 1.3, cex = cexabc, las = 1, line = 5)
mtext(expression(italic(r)[latitude-richness]), 2, cex = cexlab, line = 5)

#extra tick marks showing % of max richness
axis(1,at=log10(pcts*max(lat.corr.output$clade.richness)), labels=F,tck= .01)
text(log10(pcts*max(lat.corr.output$clade.richness)), rep(-1.12,length(pcts)), paste(pcts*100,"%",sep=""), cex = 1.75) 

dev.off()
