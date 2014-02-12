#Animation run in R window showing how the geographic pattern of richness varies thru time

source('code/supplemental_analysis_functions.r')
sim.matrix = read.csv('SENC_Master_Simulation_Matrix.csv', header = T)

# Function that shows a time lapse movie of the development of the latitudinal gradient in species richness.
# Arguments include the simID, sim.matrix, directory in which simulation files are stored, 
# temporal resolution and duration of the movie, and whether the sim files need to be unzipped or not.
# Now modified to compare two simulations at once.

lat.grad.movie.comp = function(sims, sim.matrix, time.step, time.max, unzip=F) {
  
  if(unzip) {
    sim.out1 = output.unzip(sim_dir, sims[1])
    sim.out2 = output.unzip(sim_dir, sims[2])
    all.populations1 = sim.out1$all.populations
    all.populations2 = sim.out2$all.populations
  } else { 
    all.populations1 = read.csv(paste('raw_sim_output/sim', sims[1], '_out/SENC_all.pops_sim',sims[1],'.csv',sep=''), header=T)
    all.populations2 = read.csv(paste('raw_sim_output/sim', sims[2], '_out/SENC_all.pops_sim',sims[2],'.csv',sep=''), header=T)
  }
  
  params1 = sim.matrix[sim.matrix$sim.id==sims[1],]
  params2 = sim.matrix[sim.matrix$sim.id==sims[2],]
  
  timeslices = seq(time.step, time.max, by=time.step)

  max.rich = max( c( table(all.populations1[all.populations1$time.of.extinction > 30000,'region']), 
                     table(all.populations2[all.populations2$time.of.extinction > 30000,'region'])))
  
  #reg.rich.thru.time = data.frame(time=NA, region=NA, total.rich=NA)
  par(mfrow = c(2,1), mar = c(4,4,2,1), mgp = c(2.5,1,0))
  for (t in timeslices) {
    
    all.pops1 = subset(all.populations1, time.of.origin < t & time.of.extinction > t)
    all.reg.rich1 = data.frame(table(all.pops1$region))
    names(all.reg.rich1) = c('region','total.rich')
    all.reg.rich1$region = as.numeric(as.character(all.reg.rich1$region))
    
    all.pops2 = subset(all.populations2, time.of.origin < t & time.of.extinction > t)
    all.reg.rich2 = data.frame(table(all.pops2$region))
    names(all.reg.rich2) = c('region','total.rich')
    all.reg.rich2$region = as.numeric(as.character(all.reg.rich2$region))
    
    Sys.sleep(0)
    # Panel 1: richness gradient
    plot(11 - all.reg.rich1$region, log10(all.reg.rich1$total.rich), type='b', lwd = 4, cex = .5, col = 'red',
         xlim = c(0,11), ylim=c(0, log10(max.rich)+1), xlab="Latitude",ylab = "log10 Species richness",
         main=paste(params1[1,3],"origin; K", params1[1,8], "; w =",params1[1,4],"; sigma = ",params1[1,7]))
    points(11 - all.reg.rich2$region, log10(all.reg.rich2$total.rich), type='b', lwd = 4, cex = .5, col = 'blue',
           xlim = c(0,11), ylim=c(0, log10(max.rich)+.5))
    legend('topleft',c(paste('freq =', params1[1,15],'; magn =', params1[1,16]),
                              paste('freq =', params2[1,15],'; magn =', params2[1,16])), lwd=4, col = c('red','blue'))
    axis(4, at = seq(0, ceiling(2*(log10(max.rich) + 1))/2, by = 0.5), tck = .03)
    text(10, log10(max.rich) + 1, paste("t =", t))
    
    # Panel 2: environmental optima
    regions = c(1,4,7,10)
    
    plot(1,1,type="n",xlim=c(-5,45), ylim = c(0,0.3), xlab="Thermal Optimum", ylab = "Density")
    reg.env = unique(all.populations1[all.populations1$region %in% regions, c('region','reg.env')])
    abline(v = reg.env$reg.env, lty="dotted")
    for (r in regions) {
      if (nrow(all.pops1[all.pops1$region==r,]) > 3) {
        points(density(all.pops1[all.pops1$region==r, 'env.opt']), type = 'l', col = 'red')
      }
      if (nrow(all.pops2[all.pops2$region==r,]) > 3) {
        points(density(all.pops2[all.pops2$region==r, 'env.opt']), type = 'l', col = 'blue')
      }
    }
    
    #reg.rich.thru.time = rbind(reg.rich.thru.time, cbind(time=rep(t,nrow(all.reg.rich)), all.reg.rich))
    for (i in 1:100000) {}
  }
}

# Example for plotting movie for comparing sims 4565 and 4685
lat.grad.movie.comp(sims = c(4565, 4685), sim.matrix, 500, time.max = 30000)

