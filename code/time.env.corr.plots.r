# Function for plotting 2 coefficients (corr betw richness and environment, corr betw richness
# and time) as a function of clade origin time.

 
time.env.corr.plot = function(stats.output,
                                  sim.params,
                                  min.num.data.pts = 10, 
                                  min.num.spp.per.clade = 30, 
                                  min.num.regions = 5) 
  {

  t = 30000 #max time of simulations, only used here for plotting grey rectangle
  
  if (sim.params[1,8]=='on' & sim.params[1,9]=='on') {
    K.text = 'K gradient present'
  } else if (sim.params[1,8]=='on' & sim.params[1,9]=='off') {
    K.text = 'K constant across regions'
  } else if (sim.params[1,8]=='off') {
    K.text = 'no K'
  }
    
  x = subset(stats.output, n.regions >= min.num.regions & clade.richness >= min.num.spp.per.clade)
  starting.clades = aggregate(x$clade.richness, by = list(x$sim), max)
  names(starting.clades) = c('sim','clade.richness')
  clade.origin = merge(starting.clades, x, by= c('sim','clade.richness'))
  spline.df = 4
    
  if (length(x$r.env.rich[!is.na(x$r.env.rich)]) > min.num.data.pts) {
    plot(log10(x$clade.origin.time), x$r.env.rich, xlab="",ylab="r (Env-Richness)",ylim=c(-1,1),
         main = paste("Sim ",sim.params$sim.id[1],", ",K.text,"\norigin = ",
                      sim.params[1,3],", dist.freq = ", sim.params[1,'disturb_frequency'], sep=""))
    rect(-1000,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)
    points(smooth.spline(log10(x$clade.origin.time[!is.na(x$r.env.rich)]),x$r.env.rich[!is.na(x$r.env.rich)],df=spline.df),
      type='l',col='red')
    points(log10(clade.origin$clade.origin.time), clade.origin$r.env.rich, col = 'skyblue', pch = 17, cex = 2)
  } else {
    plot(1,1,xlab="",ylab = "r (Env-Richness)",type="n",xlim=c(0,t),ylim=c(-1,1))
  }
    
  if (length(x$r.time.rich[!is.na(x$r.time.rich)]) > min.num.data.pts) {
    plot(log10(x$clade.origin.time), x$r.time.rich, xlab = "",ylab="r (Time-Richness)",ylim=c(-1,1))
    rect(-1000,.5,t+50,1.1,col=rgb(.1,.1,.1,.1),lty=0)
    points(smooth.spline(log10(x$clade.origin.time[!is.na(x$r.time.rich)]),x$r.time.rich[!is.na(x$r.time.rich)],df=spline.df),type='l',col='red')
    points(log10(clade.origin$clade.origin.time), clade.origin$r.time.rich, col = 'skyblue', pch = 17, cex = 2)
  } else { 
    plot(1,1,xlab="",ylab = "r (Time-Richness)",type="n",xlim=c(0,t),ylim=c(-1,1))
  }
  mtext("Clade origin time",1,outer=T,line=1.2)
}

