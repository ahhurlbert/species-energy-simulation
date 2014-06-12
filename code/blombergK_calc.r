# Calculate Blomberg's K for both region and environmental optimum for a given simulation
# To be run from species-energy-simulation/

BKcalc = function(simdir, sim, t=30000) {
  require(picante)
  source('Z:/git/species-energy-simulation/code/supplemental_analysis_functions.r')
  simresults = output.unzip(simdir, sim)
  all.pops = simresults$all.populations
  extant.pops = subset(all.pops, time.of.origin < t & time.of.extinction > t)
  extant.pops$spp.name = as.character(extant.pops$spp.name)
  phy = simresults$phylo.out
  extant.phy = drop.tip(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
  
  #make sure that each species only occurs once; choose the regional occurrence with best fit to env
  extant.pops$Ediff = abs(extant.pops$env.opt - extant.pops$reg.env)
  ext.pops = aggregate(extant.pops$Ediff, by = list(extant.pops$spp.name), min)
  names(ext.pops) = c('spp.name', 'Ediff')
  ext.pops2 = merge(ext.pops, extant.pops[, c('spp.name', 'region', 'env.opt', 'Ediff')], by = c('spp.name', 'Ediff'))
  
  #make sure that trait data is sorted in the same order as phy tip labels
  species = data.frame(phy$tip.label)
  names(species) = 'spp.name'
  traitdata = merge(species, ext.pops2, by = 'spp.name', sort = F)
  
  reg.BK = phylosignal(traitdata$region, extant.phy)
  env.BK = phylosignal(traitdata$env.opt, extant.phy)
  
  BKdata = rbind(reg.BK, env.BK)
  out.id = data.frame(sim = sim, var = c('region', 'env.opt'))
  outdata = cbind(out.id, BKdata)
  
  return(outdata)
}

simdirs = c('Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output',
            'Z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output',
            'raw_sim_output',
            'raw_sim_output')
sims = c(3465, 4065, 5525, 5625)

BK.output = c()
for (i in 1:4) {
  if (sims[i] == 3465) { time = 200 } else { time = 30000 }
  temp = BKcalc(simdirs[i], sims[i], t = time)
  BK.output = rbind(BK.output, temp)
}

write.csv(BK.output, 'BlombergK_4scenarios.csv', row.names=F)
