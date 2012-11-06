# Function for looking at detailed clade level patterns for 6 example clades of different age
# The 'time.buffer' is the potential +/- variation (in time units) allowed in selecting clades of a specific age


clade.exmpl.figs = function(sim.results, clade.slices=6, time.buffer=20) {
  all.pops = sim.results$all.populations
  phylo.out = sim.results$phylo.out
  max.time = max(all.pops$time.of.origin)  
  sim.params = sim.results$params.out

  
  ###############################################################################
  # Replace code here in these next 8 lines as necessary to incorporate correct
  # values for clade time of origin...
  ####################################
  #Node origin times based on phylogeny at end of simulation
  all.dist <- dist.nodes(phylo.out)
  root.dist <- all.dist[length(phylo.out$tip.label)+1, ]
  rootdist = data.frame(clade.id = names(root.dist),time.of.origin = root.dist)
    
  clade.time.slices = round(seq(max.time/clade.slices, max.time,by=max.time/clade.slices),0)

  clades = sapply(clade.time.slices, function(x) {set.seed(1); sample(rootdist$clade.id[abs(rootdist$time.of.origin - x) < time.buffer],1)})
  clds = as.numeric(as.character(clades))

  pdf(paste('clade_example_figs_sim',sim.params[1,1],'.pdf',sep=''),width=9,height=6)
  par(mfrow=c(2,3),oma = c(1,1,4,0),mar=c(4,4,1,1))
  for (i in clds) {
    cl.members = clade.members(i,phylo.out)
    cl.subset = subset(all.pops, spp.name %in% cl.members)
    cl.analysis = regional.calc(cl.subset[cl.subset<11,c('region','spp.name','time.of.origin')], phylo.out, max.time)
    attach(cl.analysis)

    if(
    plot(richness~region, xlab="Environment",ylab="Species Richness", cex=3, pch=16, col='darkblue', xlim=c(1,10))
    legend("top",legend=paste('r =',round(cor(richness,region,use="complete.obs"),2)), bty='n')
    plot(richness~MRD, xlab = "MRD", ylab="Species Richness", cex=3, pch=15, col = 'skyblue')
    legend("top",legend=paste('r =',round(cor(richness,MRD,use="complete.obs"),2)), bty='n')
    plot(richness~PSV, xlab = "PSV", ylab="Species Richness", cex=3, pch=17, col = 'darkgreen')
    legend("top",legend=paste('r =',round(cor(richness,PSV,use="complete.obs"),2)), bty='n')
    plot(richness~time.in.region, xlab = "Time in Region", ylab="Species Richness", cex=3, pch=17, col = 'yellow3')
    legend("top",legend=paste('r =',round(cor(richness,time.in.region,use="complete.obs"),2)), bty='n')
    plot(MRD~region, ylab = "MRD", xlab="Environment", cex=3, pch=18, col = 'orange', xlim=c(1,10))
    legend("top",legend=paste('r =',round(cor(region,MRD,use="complete.obs"),2)), bty='n')
    plot(PSV~region, ylab = "PSV", xlab="Environment", cex=3, pch=15, col = 'darkred', xlim=c(1,10))
    legend("top",legend=paste('r =',round(cor(region,PSV,use="complete.obs"),2)), bty='n')

    if (sim.params[8,1]==1 & sim.params[9,1]==1) {
      K.text = 'K gradient present'
    } else if (sim.params[8,1]==1 & sim.params[9,1]==2) {
      K.text = 'K constant across regions'
    } else if (sim.params[8,1]==2) {
      K.text = 'no K'
    }
    mtext(paste('Sim',sim.params[1,1],', Origin =',sim.params[3,1],', w =',sim.params[4,1],', sigma =',sim.params[7,1],
                ',disp = ',sim.params[6,1],', specn =',sim.params[5,1],',',K.text),outer=T,line=2)
    
    mtext(paste("CladeID =",i,"; Clade time of origin =", rootdist[rootdist$clade.id==i,'RD'],"; Clade richness =",length(cl.members)),3,outer=T,line=.3)
    detach(cl.analysis)
  }
  dev.off()
}