# Function for looking at detailed clade level patterns for 6 example clades of different age
# The 'time.buffer' is the potential +/- variation (in time units) allowed in selecting clades of a specific age


clade.exmpl.figs = function(sim.results, clade.slices=6, seed=1) {
  all.pops = sim.results$all.populations
  phylo.out = sim.results$phylo.out
  max.time = as.integer(max(all.pops$time.of.origin))
  sim.params = sim.results$params.out

  
  ###############################################################################
  # Replace code here in these next 8 lines as necessary to incorporate correct
  # values for clade time of origin...
  ####################################
  #Node origin times based on phylogeny at end of simulation
  all.dist <- dist.nodes(phylo.out)
  root.dist <- all.dist[length(phylo.out$tip.label)+1, ]
  rootdist = data.frame(clade.id = names(root.dist),time.of.origin = root.dist)
    
  clade.time.slices = round(seq(max.time/(clade.slices+1), max.time - max.time/(clade.slices+1),length.out=clade.slices),0)

  #Choose a clade randomly from those that are closest to the specified clade time slice,
  #specifically those that are within 10 time units of the closest clade
  clades = sapply(clade.time.slices, function(x) {set.seed(seed); 
                  sample(rootdist$clade.id[abs(rootdist$time.of.origin - x) < (min(abs(rootdist$time.of.origin - x))+10)],1)})
  clds = as.numeric(as.character(clades))

  pdf(paste('clade_example_figs_sim',sim.params[1,1],'.pdf',sep=''),width=9,height=6)
  par(mfrow=c(2,3),oma = c(1,1,4,0),mar=c(4,4,1,1))
  for (i in clds) {
    cl.members = clade.members(i,phylo.out)
    cl.subset = subset(all.pops, spp.name %in% cl.members)
    cl.analysis = regional.calc(cl.subset[c('region','spp.name','time.of.origin','reg.env','extant')], phylo.out, max.time)
    attach(cl.analysis)

    if(length(unique(richness[!is.na(richness)])) > 2) {
      plot(richness~reg.env, xlab="Environment",ylab="Species Richness", cex=3, pch=16, col='darkblue',xlim=range(sim.results$all.populations$reg.env))
      legend("top",legend=paste('r =',round(cor(richness,reg.env,use="complete.obs"),2)), bty='n')
    } else {
      plot(1,1, xlab="Environment",ylab="Species Richness", type="n")
      text(1,1,"NA")
    }
    
    if(length(unique(MRD[!is.na(MRD)])) > 2) {
      plot(richness~MRD, xlab = "MRD", ylab="Species Richness", cex=3, pch=15, col = 'skyblue')
      legend("top",legend=paste('r =',round(cor(richness,MRD,use="complete.obs"),2)), bty='n')
    } else {
      plot(1,1, xlab="MRD",ylab="Species Richness", type="n")
      text(1,1,"NA")
    }

    if(length(unique(PSV[!is.na(PSV)])) > 2) {
      plot(richness~PSV, xlab = "PSV", ylab="Species Richness", cex=3, pch=17, col = 'darkgreen')
      legend("top",legend=paste('r =',round(cor(richness,PSV,use="complete.obs"),2)), bty='n')
    } else {
      plot(1,1, xlab="PSV",ylab="Species Richness", type="n")
      text(1,1,"NA")    
    }
    
    if(length(unique(richness[!is.na(richness)])) > 2) {
      plot(richness~time.in.region, xlab = "Time in Region", ylab="Species Richness", cex=3, pch=17, col = 'yellow3')
      legend("top",legend=paste('r =',round(cor(richness,time.in.region,use="complete.obs"),2)), bty='n')
    } else {
      plot(1,1, xlab="Time in Region",ylab="Species Richness", type="n")
      text(1,1,"NA")
    }
    
    if(length(unique(MRD[!is.na(MRD)])) > 2) {
      plot(MRD~reg.env, ylab = "MRD", xlab="Environment", cex=3, pch=18, col = 'orange',xlim=range(sim.results$all.populations$reg.env))
      legend("top",legend=paste('r =',round(cor(reg.env,MRD,use="complete.obs"),2)), bty='n')
    } else {
      plot(1,1, xlab="Environment",ylab="MRD", type="n")
      text(1,1,"NA")
    }
    
    if(length(unique(MRD[!is.na(MRD)])) > 2) {
      plot(PSV~reg.env, ylab = "PSV", xlab="Environment", cex=3, pch=15, col = 'darkred',xlim=range(sim.results$all.populations$reg.env))
      legend("top",legend=paste('r =',round(cor(reg.env,PSV,use="complete.obs"),2)), bty='n')
    } else {
      plot(1,1, xlab="Environment",ylab="PSV", type="n")
      text(1,1,"NA")
    }

    if (sim.params$carry.cap=='on' & sim.params$energy.gradient=='on') {
      K.text = 'K gradient present'
    } else if (sim.params$carry.cap=='on' & sim.params$energy.gradient=='off') {
      K.text = 'K constant across regions'
    } else if (sim.params$carry.cap=='off') {
      K.text = 'no K'
    }
    mtext(paste('Sim',sim.params[1,1],', Origin =',sim.params[1,3],', w =',sim.params[1,4],', sigma =',sim.params[1,7],
                ',disp = ',sim.params[1,6],', specn =',sim.params[1,5],',',K.text),outer=T,line=2)
    
    mtext(paste("CladeID =",i,"; Clade time of origin =", rootdist[rootdist$clade.id==i,'time.of.origin'],"; Clade richness =",length(cl.members)),3,outer=T,line=.3)
    detach(cl.analysis)
  }
  dev.off()
}
