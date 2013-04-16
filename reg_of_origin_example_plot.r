
# output here is a SENC_stats_sim file that was redone (i.e., not part of any output prior
# to April 2013) with the calculation of a region of origin for each clade.
# This was originally done for sim 3325.

# The region of origin is most easily calculated using some version of this code,
# which is a modification of the clade loop in the analysis_workflow_local.r script:
regorigin = c()
for (c in (num.of.spp+1):max(sub.phylo$edge)) {
  
  #pull out list of species names belonging to each subclade
  sub.clade = clade.members(c, sub.phylo, tip.labels=T)
  subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
  
  #sub.populations is the subset of populations specific to a particular clade and timeslice
  sub.populations = subset(subset.populations, time.of.origin <= t & time.of.extinction > t)
  
  #Region of origin of clade
  reg.of.origin = sub.populations[sub.populations$time.of.sp.origin == min(sub.populations$time.of.sp.origin), 'region']
  regorigin = rbind(regorigin, cbind(sim=sim,clade.id = c, reg.of.clade.origin = reg.of.origin))
  #print(paste(sim,c,t,date(),length(sub.clade.phylo$tip.label),sep="   "));
} # end sub clade for loop


####################

#Assuming that result is merged back into the stats output, the plot is as follows

colors = colorRampPalette(c('blue','green','yellow','orange','red'))

pdf('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots/reg_of_origin_sim3325.pdf',
    height = 8, width = 5)
par(mfrow=c(2,1), mar = c(2,4,1,1), oma=c(3,1,1,1))
plot(jitter(11-output$reg.of.clade.origin), jitter(output$n.regions), pch = 16, col = colors(200)[output$colrenv],
     xlab="",ylab="Number of regions occupied",main="Energy Gradient, tropical origin, sim 3325", xaxt="n")
legend("topright",legend=c("+1","+.5","0","-.5","-1"),pch=16, col = c('blue','green','yellow','orange','red'), cex=.8)
plot(jitter(11-output$reg.of.clade.origin), log10(output$clade.richness), pch = 16, col = colors(200)[output$colrenv],
     xaxt="n",ylab="Clade richness")
mtext(c("Tropical","Temperate"),1, at = c(1,10))
mtext("Region of origin",1,outer=T)
dev.off()
