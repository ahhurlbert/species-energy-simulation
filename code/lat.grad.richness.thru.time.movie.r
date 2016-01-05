#Animation run in R window showing how the geographic pattern of richness varies thru time

source('code/supplemental_analysis_functions.r')
sim.matrix = read.csv('SENC_Master_Simulation_Matrix.csv', header = T)

# Function that shows a time lapse movie of the development of the latitudinal 
# gradient in species richness and total species richness through time.
# Arguments include the simID, sim.matrix, directory in which simulation files are stored, 
# temporal resolution and duration of the movie, and whether the sim files need to be unzipped or not.

lat.time.movie = function(sim,          #simIDs to run
                          sim_dir,      #directory where sim output is
                          movie_dir,    #directory to save movie if plot.pdf=T
                          sim.matrix,   #sim matrix with sim parameters
                          time.step,    #frequency of time lapse
                          time.max,     #final time point
                          title = F,    #put title on figure?
                          plot.pdf = F, #create pdf?
                          unzip=F) {    #do files need to be unzipped? T/F
  
                               
  if(unzip) {
    sim.out = output.unzip(sim_dir, sim[1])
    all.populations = sim.out$all.populations
  } else { 
    all.populations = read.csv(paste(sim_dir, '/sim', sim, '_out/SENC_all.pops_sim',sim,'.csv',sep=''), header=T)
  }
  
  params = sim.matrix[sim.matrix$sim.id==sim,]

  timeslices = seq(time.step, time.max, by=time.step)

  ap = unique(all.populations[, c('spp.name', 'time.of.sp.origin', 'time.of.sp.extinction')])
  rich = sapply(1:time.max, function(x) 
                length(unique(ap$spp.name[ap$time.of.sp.origin<x & ap$time.of.sp.extinction >= x])))

  if(plot.pdf) {
    pdf(paste(movie_dir, '/movie_sim',sim,'_timestep',time.step,'.pdf',sep=''), height=8, width=6, bg='white')
  }
  par(mfrow = c(2,1), mar = c(4,4,2,1), mgp = c(2.5,1,0), las = 1, cex.lab = 1.5)
  for (t in timeslices) {
    
    all.pops1 = subset(all.populations, time.of.origin < t & time.of.extinction > t)
    all.reg.rich1 = data.frame(table(all.pops1$region))
    names(all.reg.rich1) = c('region','total.rich')
    all.reg.rich1$region = as.numeric(as.character(all.reg.rich1$region))
    
    Sys.sleep(0)
    # Panel 1: richness gradient
    if(title) {
      plot.title = paste(params[1,3],"origin; K", params[1,8], "; w =",params[1,4],"; sigma = ",params[1,7])
    } else {
      plot.title = ''
    }
    
    plot(11 - all.reg.rich1$region, log10(all.reg.rich1$total.rich), type='b', lwd = 4, cex = .5, col = 'red',
         xlim = c(0,11), ylim=c(0, log10(max(rich))), ylab = "log10 Species richness", xlab = "", xaxt = "n",
         main=plot.title, cex.lab = 1.5)
    axis(4, at = seq(0, ceiling(2*(log10(max(rich))))/2, by = 0.5), tck = .03)
    text(10, log10(max(rich)), paste("t =", t))
    mtext("Latitude", 1, line = 1, cex = 1.5)
    mtext(c("Tropics","Temperate"), 1, at = c(1, 10))
    
    # Panel 2: richness through time
    plot(1:t, log10(rich[1:t]), type="l", xlim=c(0,time.max), ylim = c(0, log10(max(rich))), 
         xlab="", ylab = "log10 Species Richness", cex.lab = 1.5, xaxt = "n")
    mtext("Time", 1, line = 1, cex = 1.5)
    
    
    #reg.rich.thru.time = rbind(reg.rich.thru.time, cbind(time=rep(t,nrow(all.reg.rich)), all.reg.rich))
    for (i in 1:100000) {}
  }
  if(plot.pdf) {dev.off()}
}

#This creates a multipage pdf when plot.pdf=T. To convert this to a gif animation, install ImageMagick and GhostScript,
#both freely available. After they are installed, the following command will work.

#shell("convert -delay 50 -density 150 movie_sim3345_timestep500.pdf movie_sim3345_timestep500.gif")

# -delay specifies the delay between frames in 1/100 s
# -density specifies the dpi
# on a MacOS, use the command 'system' rather than 'shell'

# Note that the current working directory must be where the pdf file was saved.
# Also, if this command only converts the first page of the pdf to a gif, then 
# edit the delegates.xml file associated with ImageMagick, and change
# "pngalpha" to "pnmraw".

# Note that I now have to explicitly set the background color of the pdf to be "white" or else
# the conversion to a .gif leads to a black background.
