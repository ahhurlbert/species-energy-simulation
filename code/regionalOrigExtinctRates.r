

regionalRates = function(all.pops, whichRegions, timeslices) {
  glob.pops = unique(all.populations[, c('spp.name', 'region', 'time.of.origin', 'time.of.extinction')])
  
  # Calculate origination and extinction rates for a given region only for the populations
  # that originated and went extinct in that region. That is, ignore originations due to 
  # dispersal from adjacent regions and ignore extinctions of populations that originated elsewhere.
  reg.of.origin = aggregate(glob.pops$time.of.origin, by = list(glob.pops$spp.name), min)
  names(reg.of.origin) = c('spp.name', 'time.of.origin')
  pops = merge(reg.of.origin, glob.pops, by = c('spp.name', 'time.of.origin'), all.x = T)

  dataOut = c()
  for (t in timeslices) {
    t0 = max(c(0, timeslices[timeslices < t]))
    for (reg in whichRegions) {
      origin.pops = subset(pops, region == reg & 
                             time.of.origin >= t0 &
                             time.of.origin < t &
                             time.of.extinction > t)
      origins = nrow(origin.pops)
      extinct.pops = subset(pops, region == reg & 
                              time.of.origin <= t0 &
                              time.of.extinction >= t0 &
                              time.of.extinction < t)
      extincts = nrow(extinct.pops)
      time.interval = t - t0
      richness.t0 = nrow( subset(glob.pops, region == reg & time.of.origin <= t0 & time.of.extinction > t0) )
      dataOut = rbind(dataOut, c(reg, t, origins, extincts, richness.t0, time.interval))
    } #end regions loop
  } #end time loop
  dataOut = data.frame(dataOut)
  names(dataOut) = c('region', 'time', 'raw.origins', 'raw.extinctions', 'richness.t0', 'time.interval')  
}

