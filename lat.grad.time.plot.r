# This function takes simulation output and creates a plot of the latitudinal gradient in
# species richness for a specified number of points throughout the time course of the simulation.

lat.grad.time.plot = function(sim.results, numslices) {
  all.pops = sim.results$all.populations
  rich.time = sim.results$time.richness
  phylo.out = sim.results$phylo.out
  sim.params = sim.results$end.params
  
  #Eliminate results from boundary bins (0, 11)
  rich.time2 = rich.time[!rich.time$region %in% c(0,11),]
  max.time = max(rich.time$time)

  slices = round(seq(max.time/numslices, max.time,by=max.time/numslices),0)
  cols = rainbow(numslices+1)

  pdf(paste('lat_grad_thru_time_sim',sim,'.pdf',sep=''),height=6,width=8)
  par(mfrow=c(1,1),oma=c(0,0,2.5,0),mar=c(4,4,1,1))
  plot(c(1,10),c(0,max(rich.time2$spp.rich)),type="n",xlab="Latitude",ylab="Species richness")
  sapply(slices,function(x) points(11 - rich.time2$region[rich.time2$time==x], rich.time2$spp.rich[rich.time2$time==x],type='l',col=cols[which(slices==x)]))
  legend('topright',legend=slices,lty='solid',col=cols[2:(numslices+1)])
  if (sim.params[8,1]==1 & sim.params[9,1]==1) {
    K.text = 'K gradient present'
  } else if (sim.params[8,1]==1 & sim.params[9,1]==2) {
    K.text = 'K constant across regions'
  } else if (sim.params[8,1]==2) {
    K.text = 'no K'
  }
  mtext(paste('Sim',sim.params[1,1],', Origin =',sim.params[3,1],', w =',sim.params[4,1],', sigma =',sim.params[7,1],
              ',\ndisp = ',sim.params[6,1],', specn =',sim.params[5,1],',',K.text),outer=T)
  dev.off()
}

