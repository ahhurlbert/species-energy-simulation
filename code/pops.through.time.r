


# calculate species richness and total cumulative number of population extinctions and initiations
# this is done within each bin at a given timeslice

pops.through.time = function(all.populations, timeslice) {
  sub.populations = subset(all.populations, time.of.origin <= timeslice)
  num.all.lineages = data.frame(table(sub.populations$region))
  extinct.lineages = aggregate(sub.populations$time.of.extinction,by=list(sub.populations$region), function(x) sum(x<timeslice))
  time.col = aggregate(sub.populations$time.of.origin,by=list(sub.populations$region), min)
  output = data.frame(region = extinct.lineages$Group.1, extinct.pops=extinct.lineages$x, 
                      total.pops=num.all.lineages$Freq, time.in.region.pops=timeslice-time.col$x,extant.pops = num.all.lineages$Freq - extinct.lineages$x)
  return(output)   
}

