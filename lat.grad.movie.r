github_dir = 'C:/Documents and Settings/Hurlbert/species-energy-simulation'
source(paste(github_dir,'/unzipping_files.r',sep=''))

sim_dir = 'c:/sencoutput/senc.out.130115' #directory with simulation output

#sim_matrix = 'c:/sencoutput/senc.out.130115/sim.matrix.output_2013-01-16.csv'
sim.matrix = read.csv('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/sim.matrix.output_2012-12-14.csv',header=T)


# Function that shows a time lapse movie of the development of the latitudinal gradient in species richness.
# Arguments include the simID, sim.matrix, directory in which simulation files are stored, 
# temporal resolution and duration of the movie, and whether the sim files need to be unzipped or not.
lat.grad.movie = function(sim, sim.matrix, sim_dir, time.step, time.max, unzip=T) {
  
  if(unzip) {
    sim.out = output.unzip(sim_dir, sim)
    all.populations = sim.out$all.populations
  } else { 
    all.populations = read.csv(paste(sim_dir,'/SENC_all.pops_sim',sim,'.csv',sep=''), header=T)
  }
  
  params = sim.matrix[sim.matrix$sim.id==sim,]
  
  timeslices = seq(time.step, time.max, by=time.step)

  max.rich = max(table(all.populations[all.populations$time.of.extinction > 30000,'region']))
  
  reg.rich.thru.time = data.frame(time=NA, region=NA, total.rich=NA)
  par(mfrow = c(1,1))
  for (t in timeslices) {
    
    all.pops = subset(all.populations, time.of.origin < t & time.of.extinction > t)
    all.reg.rich = data.frame(table(all.pops$region))
    names(all.reg.rich) = c('region','total.rich')
    all.reg.rich$region = as.numeric(as.character(all.reg.rich$region))

    plot(11 - all.reg.rich$region, log10(all.reg.rich$total.rich), type='b', lwd = 4, cex = .5, col = 'red',
         xlim = c(0,11), ylim=c(0, log10(max.rich)+.5), xlab="Latitude",ylab = "log10 Species richness",
         main=paste(params[1,3],"origin; K", params[1,8], "; w =",params[1,4],"; sigma = ",params[1,7]))
    text(10, log10(max.rich)+.5, paste("t =", t))
    
    reg.rich.thru.time = rbind(reg.rich.thru.time, cbind(time=rep(t,nrow(all.reg.rich)), all.reg.rich))
    for (i in 1:1000000) {}
  }
}

# Example for plotting movie for sim 2635
lat.grad.movie(2635, sim_dir, 10)

