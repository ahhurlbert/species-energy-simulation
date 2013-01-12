# Function for plotting 3 coefficients (corr betw richness and environment, corr betw richness
# and time, gamma) as a function of clade origin time.

# NOTE: May not be worth doing this over early time slices. Perhaps just for last time point,
#       in which case time loop can be commented out and set t = max(stats.output$time)
 
clade.origin.corr.plot.simple = function(stats.output,
                                  sim.params,
                                  min.num.data.pts = 10, 
                                  min.num.spp.per.clade = 30, 
                                  min.num.regions = 5, 
                                  output.dir) {
  timeslices = unique(stats.output$time[!is.na(stats.output$time)])  

  
  for (t in timeslices) {
    #t = max(stats.output$time)                 #comment this out if looping over multiple time slices
    x = subset(stats.output, time==t & n.regions >= min.num.regions & clade.richness >= min.num.spp.per.clade, select = 2:ncol(stats.output))
    spline.df = 4
    
    if (length(x$r.env.rich[!is.na(x$r.env.rich)]) > min.num.data.pts) {
      plot(x$clade.origin.time, x$r.env.rich, xlab="",ylab="Environment-Richness correlation",ylim=c(-1,1),main=paste('w =',sim.params[1,4],', sigma =',sim.params[1,7],sep=""))
      points(smooth.spline(x$clade.origin.time[!is.na(x$r.env.rich)],x$r.env.rich[!is.na(x$r.env.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1000,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.time.rich[!is.na(x$r.time.rich)]) > min.num.data.pts) {
      plot(x$clade.origin.time, x$r.time.rich, xlab = "",ylab="Time in region-Richness correlation",ylim=c(-1,1),main=paste('w =',sim.params[1,4],', sigma =',sim.params[1,7],sep=""))
      points(smooth.spline(x$clade.origin.time[!is.na(x$r.time.rich)],x$r.time.rich[!is.na(x$r.time.rich)],df=spline.df),type='l',col='red')
    } else { 
      plot(1,1,xlab="",ylab = "Time in Region-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1000,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    mtext("Clade origin time",1,outer=T,line=2)
    if (sim.params[1,8]=='on' & sim.params[1,9]=='on') {
      K.text = 'K gradient present'
    } else if (sim.params[1,8]=='on' & sim.params[1,9]=='off') {
      K.text = 'K constant across regions'
    } else if (sim.params[1,8]=='off') {
      K.text = 'no K'
    }
    mtext(paste('Origin =',sim.params[1,3],', ndisp = ',sim.params[1,6],', specn =',sim.params[1,5],',',K.text),outer=T)
  
    if (length(x$gamma.stat[!is.na(x$gamma.stat)]) > min.num.data.pts) {
      plot(x$clade.origin.time, x$gamma.stat, xlab = "", ylab="Gamma",main=paste('w =',sim.params[1,4],', sigma =',sim.params[1,7],sep=""))
      points(smooth.spline(x$clade.origin.time[!is.na(x$gamma.stat)],x$gamma.stat[!is.na(x$gamma.stat)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Gamma",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1000,-1.645,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0) # -1.645 is the critical value for rejecting constant rates

  }
 
  
}

