# abundance_thru_time.r
#
# function for plotting simulated abundances within a particular region through time
env.fit.fun = function (w.param,E,spp.opt) { 
  exp(-(spp.opt - E)^2/(w.param^2)) 
}

realPopSize = function(all.pops, time, suppress.warning = F) {
  pops.extant = subset(all.pops, time.of.extinction > time & time.of.origin < time, 
                       select = c('spp.name', 'region', 'env.opt', 'reg.env', 'time.of.sp.origin',
                                  'time.of.origin', 'time.of.extinction', 'carry.cap', 'rawPopSize'))
  if (nrow(pops.extant) == 0 & suppress.warning == F) { 
    print('No populations present at this time') ; return(NULL)
  } else if (nrow(pops.extant) == 0 & suppress.warning == T) { 
    return(NULL) 
  }
  
  #The pop.multiplier assures that the sum of species abundances will not exceed carrying capacity
  pops.extant$pop.multiplier = pops.extant$carry.cap / sum(pops.extant$rawPopSize)
  pops.extant$pop.multiplier[pops.extant$pop.multiplier > 1] = 1
  pops.extant$realPopSize = round( pops.extant$rawPopSize * pops.extant$pop.multiplier , digits = 0)
  pops.extant$time = time
  return(pops.extant)
}

# Within a specified region, calculate individual species abundances at each point in time
# over a particular time series
abundanceThruTime = function(simData, Region, startTime, endTime, timeStep, dataFrameOut = T, dots = T, richnessPlot = T) {
  params = simData$params.out
  all.pops = simData$all.populations
  time.rich = simData$time.richness
  
  regPops = subset(all.pops, region == Region)
  
  regPops$rawPopSize = env.fit.fun(params$w, regPops$reg.env, regPops$env.opt) * regPops$carry.cap
  
  times = seq(startTime, endTime, by = timeStep)
  
  #abundances
  popsAbundThruTime = c()
  for (t in times) {
    popsAbund = realPopSize(regPops, t, suppress.warning = T)
    popsAbundThruTime = rbind(popsAbundThruTime, popsAbund)
  }
  
  # Abundance trajectories
  par(mar = c(4,4,4,4), las = 1)
  plot(range(times), range(log10(popsAbundThruTime$realPopSize)), type = "n", xlab = "Time", 
       ylab = expression(paste(plain(log)[10]," Abundance")), main = paste("Sim", params$sim.id))
  uniq.spp = unique(popsAbundThruTime$spp.name)
  sapply(uniq.spp, function(x) { temp = subset(popsAbundThruTime, spp.name == x); 
                                 points(temp$time, log10(temp$realPopSize), type = 'l')})
  
   # Plot dots indicating origins and extinctions if dots == T
  if (dots) {
    # Place small red dots at the time of last occurrence for each species lasting longer than one timestep
    # and small blue dots at the first time of occurrence
    time.extinct = aggregate(popsAbundThruTime$time, by = list(popsAbundThruTime$spp.name), max)
    names(time.extinct) = c('spp.name', 'time')
    time.origin = aggregate(popsAbundThruTime$time, by = list(popsAbundThruTime$spp.name), min)
    names(time.origin) = c('spp.name', 'time')
    durableSpp = time.extinct$spp.name[time.extinct$time - time.origin$time > timeStep]
    sapply(durableSpp, function(x) { 
      preExtinctAbund = popsAbundThruTime$realPopSize[popsAbundThruTime$spp.name == x & 
                                                     popsAbundThruTime$time == time.extinct$time[time.extinct$spp.name == x]]
      points(time.extinct$time[time.extinct$spp.name == x], log10(preExtinctAbund), pch = 16, cex = .5, col = 'red') 
      originalAbund = popsAbundThruTime$realPopSize[popsAbundThruTime$spp.name == x & 
                                                   popsAbundThruTime$time == time.origin$time[time.origin$spp.name == x]]
      points(time.origin$time[time.origin$spp.name == x], log10(originalAbund), pch = 16, cex = .5, col = 'blue') })
  }
  # Richness trajectory
  if (richnessPlot) {
    regTimeRich = subset(time.rich, region == Region & time %in% times)
    par(new = T)
    plot(regTimeRich$time, regTimeRich$spp.rich, xaxt = "n", yaxt = "n", 
         xlab = "", ylab = "", type = 'l', lwd = 5, col = 'gray90')
    axis(4)
    mtext("Richness", 4, las = 0, line = 3)
  }
  if (dataFrameOut) {
    return(popsAbundThruTime)  
  }
}