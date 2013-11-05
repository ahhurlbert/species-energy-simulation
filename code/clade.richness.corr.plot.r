# Function for plotting the 6 correlation coefficients (e.g., betw richness and MRD, richness
# and time, MRD and environment,etc) as a function of clade origin time.

# NOTE: May not be worth doing this over early time slices. Perhaps just for last time point,
#       in which case time loop can be commented out and set t = max(stats.output$time)
 
clade.richness.corr.plot = function(stats.output,
                                  sim.params,
                                  min.num.data.pts = 10, 
                                  min.num.spp.per.clade = 30, 
                                  min.num.regions = 5, 
                                  output.dir) {
  timeslices = unique(stats.output$time[!is.na(stats.output$time)])  

  #Plotting
  pdf(paste(output.dir,'/corrs_vs_claderichness_sim',stats.output$sim[2],'.pdf',sep=''),height=8,width=9)
  par(mfrow=c(3,3),oma=c(5,1,4,0),mar=c(2,5,2,1))

  for (t in timeslices) {
    #t = max(stats.output$time)                 #comment this out if looping over multiple time slices
    x = subset(stats.output, time==t & n.regions >= min.num.regions & clade.richness >= min.num.spp.per.clade, select = 2:ncol(stats.output))
    spline.df = 4
    
    if (length(x$r.env.rich[!is.na(x$r.env.rich)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.env.rich, xlab="",ylab="Environment-Richness correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.env.rich)]),x$r.env.rich[!is.na(x$r.env.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,.5,log10(max(x$clade.richness))+.5,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.MRD.rich[!is.na(x$r.MRD.rich)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.MRD.rich, xlab="",ylab="Mean Root Distance-Richness correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.MRD.rich)]),x$r.MRD.rich[!is.na(x$r.MRD.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Mean Root Distance-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,-1.1,log10(max(x$clade.richness))+.5,-.5,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.PSV.rich[!is.na(x$r.PSV.rich)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.PSV.rich, xlab = "", ylab="PSV-Richness correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.PSV.rich)]),x$r.PSV.rich[!is.na(x$r.PSV.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "PSV-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,.5,log10(max(x$clade.richness))+.5,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.time.rich[!is.na(x$r.time.rich)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.time.rich, xlab = "",ylab="Time in region-Richness correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.time.rich)]),x$r.time.rich[!is.na(x$r.time.rich)],df=spline.df),type='l',col='red')
    } else { 
      plot(1,1,xlab="",ylab = "Time in Region-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,.5,log10(max(x$clade.richness))+.5,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.env.MRD[!is.na(x$r.env.MRD)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.env.MRD, xlab="", ylab="Environment-MRD correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.env.MRD)]),x$r.env.MRD[!is.na(x$r.env.MRD)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-MRD correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,-1.1,log10(max(x$clade.richness))+.5,-.5,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.env.PSV[!is.na(x$r.env.PSV)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.env.PSV, xlab="", ylab="Environment-PSV correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.env.PSV)]),x$r.env.PSV[!is.na(x$r.env.PSV)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-PSV correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,.5,log10(max(x$clade.richness))+.5,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.rich.ext[!is.na(x$r.rich.ext)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.rich.ext, xlab = "", ylab="Extinction Rate-Richness correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.rich.ext)]),x$r.rich.ext[!is.na(x$r.rich.ext)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Extinction Rate-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    #rect(-1000,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(x$r.ext.reg[!is.na(x$r.ext.reg)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$r.ext.reg, xlab = "", ylab="Extinction Rate-Region correlation",ylim=c(-1,1))
      points(smooth.spline(log10(x$clade.richness[!is.na(x$r.ext.reg)]),x$r.ext.reg[!is.na(x$r.ext.reg)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Extinction Rate-Region correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    #rect(-1000,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)
    
    if (length(x$gamma.stat[!is.na(x$gamma.stat)]) > min.num.data.pts) {
      plot(log10(x$clade.richness), x$gamma.stat, xlab = "", ylab="Gamma")
      points(smooth.spline(log10(x$clade.richness[!is.na(x$gamma.stat)]),x$gamma.stat[!is.na(x$gamma.stat)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Gamma",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-1,-1.645,log10(max(x$clade.richness))+.5,1.1,col=rgb(.1,.1,.1,.1),lty=0) # -1.645 is the critical value for rejecting constant rates

    mtext("Clade Richness",1,outer=T,line=2)
    if (sim.params[1,8]=='on' & sim.params[1,9]=='on') {
      K.text = 'K gradient present'
    } else if (sim.params[1,8]=='on' & sim.params[1,9]=='off') {
      K.text = 'K constant across regions'
    } else if (sim.params[1,8]=='off') {
      K.text = 'no K'
    }
    mtext(paste('Sim',sim.params[1,1],', Origin =',sim.params[1,3],', w =',sim.params[1,4],', sigma =',
                sim.params[1,7],', dist.freq =',sim.params[1,'disturb_frequency'],', dist.intens.temp =',
                sim.params[1,'temperate_disturb_intensity'],
                ',\ndisp = ',sim.params[1,6],', specn =',sim.params[1,5],',',K.text,', time =',t),outer=T)
    
  }
  dev.off()
}

