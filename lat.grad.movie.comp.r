Allen = 0;

if (Allen == 1) {
  github_dir = 'C:/Documents and Settings/Hurlbert/species-energy-simulation'
  sim_dir = 'c:/sencoutput/sims.out.130204' #directory with simulation output
}

if (Allen == 0) {
  github_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation'
  sim_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204' #directory with simulation output
}

source(paste(github_dir,'/unzipping_files.r',sep=''))
sim.matrix = read.csv(paste(github_dir,'/SENC_Master_Simulation_Matrix.csv',sep=''), header=T)

# Function that shows a time lapse movie of the development of the latitudinal gradient in species richness.
# Arguments include the simID, sim.matrix, directory in which simulation files are stored, 
# temporal resolution and duration of the movie, and whether the sim files need to be unzipped or not.
# Now modified to compare two simulations at once.

lat.grad.movie.comp = function(sims, sim.matrix, sim_dir, time.step, time.max, unzip=T) {
  
  if(unzip) {
    sim.out1 = output.unzip(sim_dir, sims[1])
    sim.out2 = output.unzip(sim_dir, sims[2])
    all.populations1 = sim.out1$all.populations
    all.populations2 = sim.out2$all.populations
  } else { 
    all.populations1 = read.csv(paste(sim_dir,'/SENC_all.pops_sim',sims[1],'.csv',sep=''), header=T)
    all.populations2 = read.csv(paste(sim_dir,'/SENC_all.pops_sim',sims[2],'.csv',sep=''), header=T)
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
         xlim = c(0,11), ylim=c(0, log10(max.rich)+.5), xlab="Latitude",ylab = "log10 Species richness",
         main=paste(params1[1,3],"origin; K", params1[1,8], "; w =",params1[1,4],"; sigma = ",params1[1,7]))
    points(11 - all.reg.rich2$region, log10(all.reg.rich2$total.rich), type='b', lwd = 4, cex = .5, col = 'blue',
           xlim = c(0,11), ylim=c(0, log10(max.rich)+.5))
    legend('topleft',c(paste('freq =', params1[1,15],'; magn =', params1[1,16]),
                              paste('freq =', params2[1,15],'; magn =', params2[1,16])), lwd=4, col = c('red','blue'))
    text(10, log10(max.rich)+.5, paste("t =", t))
    
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

# Example for plotting movie for comparing sims 3035 and 3095
lat.grad.movie.comp(sims = c(3035, 3095), sim.matrix, sim_dir, 50,time.max=30000)

