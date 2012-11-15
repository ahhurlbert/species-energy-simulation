# Function for plotting the 6 correlation coefficients (e.g., betw richness and MRD, richness
# and time, MRD and environment,etc) as a function of clade origin time.

# NOTE: May not be worth doing this over early time slices. Perhaps just for last time point,
#       in which case time loop can be commented out and set t = max(stats.output$time)

clade.origin.corr.plot = function(stats.output,sim.params,min.num.data.pts = 10) {
  timeslices = unique(stats.output$time)  

  #Plotting
  pdf(paste('corrs_vs_cladeage_sim',stats.output$sim[1],'.pdf',sep=''),height=6,width=9)
  par(mfrow=c(2,4),oma=c(5,1,4,0),mar=c(2,5,2,2))

  for (t in timeslices) {
    #t = max(stats.output$time)                 #comment this out if looping over multiple time slices
    colnames(stats.output)[1] = 'sim.id'
    stats2x = subset(stats.output, time==t)
    attach(stats2x)
    spline.df = 4
    
    if (length(r.env.rich[!is.na(r.env.rich)]) > min.num.data.pts) {
      plot(clade.origin.time, r.env.rich, xlab="",ylab="Environment-Richness correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.env.rich)],r.env.rich[!is.na(r.env.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-50,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(r.MRD.rich[!is.na(r.MRD.rich)]) > min.num.data.pts) {
      plot(clade.origin.time, r.MRD.rich, xlab="",ylab="Mean Root Distance-Richness correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.MRD.rich)],r.MRD.rich[!is.na(r.MRD.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Mean Root Distance-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-50,-1.1,t+50,-.5,col=rgb(.1,.1,.1,.1),lty=0)

    
    if (length(r.PSV.rich[!is.na(r.PSV.rich)]) > min.num.data.pts) {
      plot(clade.origin.time, r.PSV.rich, xlab = "", ylab="PSV-Richness correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.PSV.rich)],r.PSV.rich[!is.na(r.PSV.rich)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "PSV-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-50,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(r.PSV.rich[!is.na(r.PSV.rich)]) > min.num.data.pts) {
      plot(clade.origin.time, r.rich.ext, xlab = "", ylab="Extinction Rate-Richness correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.rich.ext)],r.rich.ext[!is.na(r.rich.ext)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Extinction Rate-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    #rect(-50,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(r.time.rich[!is.na(r.time.rich)]) > min.num.data.pts) {
      plot(clade.origin.time, r.time.rich, xlab = "",ylab="Time in region-Richness correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.time.rich)],r.time.rich[!is.na(r.time.rich)],df=spline.df),type='l',col='red')
    } else { 
      plot(1,1,xlab="",ylab = "Time in Region-Richness correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-50,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    
    if (length(r.env.MRD[!is.na(r.env.MRD)]) > min.num.data.pts) {
      plot(clade.origin.time, r.env.MRD, xlab="", ylab="Environment-MRD correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.env.MRD)],r.env.MRD[!is.na(r.env.MRD)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-MRD correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-50,-1.1,t+50,-.5,col=rgb(.1,.1,.1,.1),lty=0)

    
    if (length(r.env.PSV[!is.na(r.env.PSV)]) > min.num.data.pts) {
      plot(clade.origin.time, r.env.PSV, xlab="", ylab="Environment-PSV correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.env.PSV)],r.env.PSV[!is.na(r.env.PSV)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Environment-PSV correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    rect(-50,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)

    if (length(r.ext.reg[!is.na(r.ext.reg)]) > min.num.data.pts) {
      plot(clade.origin.time, r.ext.reg, xlab = "", ylab="Extinction Rate-Region correlation")
      points(smooth.spline(clade.origin.time[!is.na(r.ext.reg)],r.ext.reg[!is.na(r.ext.reg)],df=spline.df),type='l',col='red')
    } else {
      plot(1,1,xlab="",ylab = "Extinction Rate-Region correlation",type="n",xlim=c(0,t),ylim=c(-1,1))
    }
    #rect(-50,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)
    
    mtext("Clade origin time",1,outer=T,line=2)
    if (sim.params[1,8]=='on' & sim.params[1,9]=='on') {
      K.text = 'K gradient present'
    } else if (sim.params[1,8]=='on' & sim.params[1,9]=='off') {
      K.text = 'K constant across regions'
    } else if (sim.params[1,8]=='off') {
      K.text = 'no K'
    }
    mtext(paste('Sim',sim.params[1,1],', Origin =',sim.params[1,3],', w =',sim.params[1,4],', sigma =',sim.params[1,7],
                ',\ndisp = ',sim.params[1,6],', specn =',sim.params[1,5],',',K.text,', time =',t),outer=T)
    
    detach(stats2x)
  }
  dev.off()
}

