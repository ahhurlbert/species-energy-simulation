source('z:/git/species-energy-simulation/code/supplemental_analysis_functions.r')
source('z:/git/species-energy-simulation/code/abundance_thru_time.r')
sim4065 = output.unzip('z:/sencoutput/hurlbert_and_stegen_2014/raw_sim_output', 4065)
all.pops = sim4065$all.populations
extant.pops = subset(all.pops, extant==1)

#Range size distributions (trivial since few species occur outside one bin)
range.sizes = data.frame(table(extant.pops$spp.name))
names(range.sizes) = c('spp.name', 'rangesize')

#Age distributions
ages.mean = aggregate(extant.pops$time.of.sp.origin, by = list(extant.pops$region), function(x) 
                  (mean(max(extant.pops$time.of.extinction) - x)))
ages.sd = aggregate(extant.pops$time.of.sp.origin, by = list(extant.pops$region), function(x) 
  (var(max(extant.pops$time.of.extinction) - x))^0.5)

#Abundance distributions

meanAges = function(path, sim) {
  simresults = output.unzip(path, sim)
  all.pops = simresults$all.populations
  extant.pops = subset(all.pops, extant == 1)
  ages.mean = aggregate(extant.pops$time.of.sp.origin, by = list(extant.pops$region), function(x)
                        mean(max(simresults$time.richness$time) - x))
  return(ages.mean)
}