Allen = 1;

if (Allen == 1) {
  github_dir = 'C:/Documents and Settings/Hurlbert/species-energy-simulation'
  sim_dir = 'c:/sencoutput/' #directory with simulation output
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
lat.grad.movie = function(sim, sim.matrix, sim_dir, time.step, time.max, unzip=F, plot.pdf=F) {
  
  if(unzip) {
    sim.out = output.unzip(sim_dir, sim)
    all.populations = sim.out$all.populations
  } else { 
    all.populations = read.csv(paste(sim_dir,'/SENC_all.pops_sim',sim,'.csv',sep=''), header=T)
  }
  
  params = sim.matrix[sim.matrix$sim.id==sim,]
  
  timeslices = seq(time.step, time.max, by=time.step)

  max.rich = max(table(all.populations[all.populations$time.of.extinction > 30000,'region']))

  # Regions for highlighting
  regions = c(1,4,7,10)
  reg.cols = c('lightslateblue', 'lightgreen', 'orange', 'red3')
    
  reg.rich.thru.time = data.frame(time=NA, region=NA, total.rich=NA)
  if(plot.pdf) {
    pdf(paste(sim_dir,'/movies/movie_sim',sim,'_timestep',time.step,'.pdf',sep=''), height=8, width=8, bg='white')
  }
  par(mfrow = c(2,1), mar = c(4,4,1,1), mgp = c(2.5,1,0))
  
  for (t in timeslices) {
    
    all.pops = subset(all.populations, time.of.origin < t & time.of.extinction > t)
    all.reg.rich = data.frame(table(all.pops$region))
    names(all.reg.rich) = c('region','total.rich')
    all.reg.rich$region = as.numeric(as.character(all.reg.rich$region))
    Sys.sleep(0)
    plot(all.reg.rich$region, log10(all.reg.rich$total.rich), type='b', lwd = 4, cex = .5, xlab = "",
         xlim = c(0,11), ylim=c(0, log10(max.rich)+.5), ylab = expression("log"[10]*" Species richness"), xaxt = "n",
         main=paste(params[1,3],"origin; K", params[1,8], "; w =",params[1,4],"; sigma = ",params[1,7]), las = 1)
    points(all.reg.rich$region[all.reg.rich$region %in% regions], 
           log10(all.reg.rich$total.rich[all.reg.rich$region %in% regions]), 
           pch = 16, cex = 2, col = reg.cols)
    text(10, log10(max.rich)+.3, paste("t =", t))
    mtext(c("Temperate","Tropics"),1,at=c(1.5,9.5),line=0.5, cex = 1.5)
    
    reg.rich.thru.time = rbind(reg.rich.thru.time, cbind(time=rep(t,nrow(all.reg.rich)), all.reg.rich))
    
    # Panel 2: environmental optima
    plot(1,1,type="n",xlim=c(0,40), ylim = c(0,0.4), xlab="Species' Thermal Optimum (°C)", ylab = "Density", las = 1)
    reg.env = unique(all.populations[all.populations$region %in% regions, c('region','reg.env')])
    abline(v = reg.env$reg.env, lty="dotted", col = reg.cols, lwd = 2)
    for (r in regions) {
      if (nrow(all.pops[all.pops$region==r,]) > 3) {
        points(density(all.pops[all.pops$region==r, 'env.opt']), type = 'l', col = reg.cols[r==regions], lwd = 4)
      }
    }
    
    for (i in 1:100000) {}
  }
  if(plot.pdf) {dev.off()}
}

# Example for plotting movie for sim 3345
lat.grad.movie(3345, sim.matrix, sim_dir, time.step=500, time.max=30000, plot.pdf=F)

#This creates a multipage pdf when plot.pdf=T. To convert this to a gif animation, install ImageMagick and GhostScript,
#both freely available. After they are installed, the following command will work.

shell("convert -delay 50 -density 150 movie_sim3345_timestep500.pdf movie_sim3345_timestep500.gif")

# -delay specifies the delay between frames in 1/100 s
# -density specifies the dpi
# on a MacOS, use the command 'system' rather than 'shell'

# Note that the current working directory must be where the pdf file was saved.
# Also, if this command only converts the first page of the pdf to a gif, then 
# edit the delegates.xml file associated with ImageMagick, and change
# "pngalpha" to "pnmraw".

# Note that I now have to explicitly set the background color of the pdf to be "white" or else
# the conversion to a .gif leads to a black background.
