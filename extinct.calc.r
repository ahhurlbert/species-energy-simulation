


#calculate extinction rate
# This is in units of extinctions per population per time, calculated as the 
# total number of extinctions in a region divided by the total number of populations 
# that ever appeared (whether through speciation or dispersal) divided by the time
# since the first colonization of that region.
# The argument should be an all.populations dataframe subsetted down to the appropriate
# timeslice (i.e., without species with time.of.origin > t)
extinct.calc = function(all.populations, timeslice) {
  sub.populations = subset(all.populations, time.of.origin < timeslice)
  num.all.lineages = data.frame(table(sub.populations$region))
  extinct.lineages = aggregate(sub.populations$time.of.extinction,by=list(sub.populations$region), function(x) sum(x<timeslice))
  time.col = aggregate(sub.populations$time.of.origin,by=list(sub.populations$region), min)
  output = data.frame(region = extinct.lineages$Group.1, extinct.pops=extinct.lineages$x, 
                      total.pops=num.all.lineages$Freq, time.in.region.pops=timeslice-time.col$x)
  output$extinction.rate = output$extinct.pops/output$total.pops/output$time.in.region.pops
  return(output)   
}

