# Function that takes a simulation results object and a point in time (measured in time since the root)
# and slices the simulation phylogeny to only depict species extant at that time point.
# Also returned is the equivalent all.populations dataframe at that time point.

time.slice.phylo = function(simresults, t, memLimit = 4095) {
  require(ape)
  require(paleotree)
  
  #increase memory limit (for timeSliceTree below)
  memory.limit(size = memLimit)
  
  phylo.out = simresults$phylo.out
  all.populations = simresults$all.populations
  max.time.actual = max(simresults$time.richness$time)
  
  # vector of species in existence at time t
  sub.species = as.character(unique(subset(all.populations, time.of.origin <= t & time.of.extinction > t, select = 'spp.name'))[,1])
  
  # Some species may be extant globally (extant==1) but in our boundary regions (0 or 11) only;
  # we need to eliminate species that are not extant within regions 1-10 (which is all that is
  # reflected in the all.populations dataframe)
  time.slice.populations = all.populations
  time.slice.populations$extant = 0
  time.slice.populations$extant[time.slice.populations$time.of.origin <= t & time.slice.populations$time.of.extinction > t] = 1
  extant.ornot = aggregate(time.slice.populations$extant, by = list(time.slice.populations$spp.name), sum)
  extinct.species = as.character(extant.ornot[extant.ornot$x == 0, 'Group.1'])
  
  sub.species2 = sub.species[!sub.species %in% extinct.species] # 'true' list of extant species
  tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]) # species to drop
  
  # Drop extinct species out of the phylogeny and tidy up
  sub.phylo = drop.tip(phylo.out, tips.to.drop)
  temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo), Ntip(sub.phylo) + 1])
  most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo), Ntip(sub.phylo) + 1])))]
  extinct.time.most.recent = unique(all.populations$time.of.sp.extinction[all.populations$spp.name == most.recent.spp])
  sub.phylo$root.time = temp.root.time + max(c(0, max.time.actual - extinct.time.most.recent))
  sub.phylo = collapse.singles(timeSliceTree(sub.phylo, sliceTime = (max.time.actual - t), plot = F, drop.extinct = T))
  return(list(slicedphylo = sub.phylo,
              slicedpops = subset(time.slice.populations, spp.name %in% sub.species2)))
}