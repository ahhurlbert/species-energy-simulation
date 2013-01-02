# Evaluating relationship between simulation output correlations and simulation parameters

file_dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots'

#sim.rs = compile.firstlines(file_dir , "SENC_Stats_sim")          #takes ~12 minutes on Hurlbert office machine
#write.csv(sim.rs,paste(file_dir,'/allsims_bigclade_stats_output.csv',sep=''), row.names=F)

sim.rs = read.csv(paste(file_dir,'/allsims_bigclade_stats_output.csv',sep=''), header=T)
sim.matrix = read.csv(paste(file_dir,'/sim.matrix.output_2012-12-20.csv',sep=''), header=T)
# warning: sim.matrix.output_2012-12-20.csv currently messed up (all text entries changed to '1')
# I replaced first 11 columns with older sim.matrix.output file

# Set alpha and beta to explore variation in other axes of parameter space; 
Alpha = 1e-6
Beta = 1e-4

# exclude case with energy gradient but no carrying capacity
sim.matrix.sm = subset(sim.matrix, alpha == Alpha & beta == Beta & !(carry.cap=="off" & energy.gradient=="on"), 
                       select = c('sim.id','reg.of.origin','w','sigma_E','carry.cap','energy.gradient'))

sim.big = merge(sim.matrix.sm, sim.rs, by.x="sim.id", by.y="sim", all.x=T)
sim.big = sim.big[sim.big$sim.id!=0,]

# Within each combination of w and sigma, space out simulations along x-axis based
# on carry.cap/energy.gradient and along y-axis based on reg.of.origin
x.offset = .5
y.offset = .25

f1 = function(x) {
  if(sim.big[x,'carry.cap']=="on" & sim.big[x,'energy.gradient']=="on") {
     y = sim.big[x,'w'] - x.offset
  }
  if(sim.big[x,'carry.cap']=="off" & sim.big[x,'energy.gradient']=="off") {
    y = sim.big[x,'w'] + x.offset
  }
  if (sim.big[x,'carry.cap']=="on" & sim.big[x,'energy.gradient']=="off") {
    y = sim.big[x,'w']
  }
  return(y)
}
sim.big$w.K = sapply(1:nrow(sim.big), f1)

f2 = function(x) {
  if(sim.big[x,'reg.of.origin']=='tropical') {
    y = sim.big[x,'sigma_E'] + .25
  }
  if(sim.big[x,'reg.of.origin']=='temperate') {
    y = sim.big[x,'sigma_E'] - .25
  }
  return(y)
}
sim.big$sigma.reg = sapply(1:nrow(sim.big), f2)

sim.big$symbol = 16
sim.big[sim.big$reg.of.origin=='temperate','symbol'] = 17

# List of independent variables to plot (using color)
yvars = c('r.time.rich','r.env.rich','r.MRD.rich','r.PSV.rich','r.env.MRD',
          'r.env.PSV','r.ext.reg','r.rich.ext')


# Plots are on sigma_E versus w space, with color reflecting yvars (list above)
# TO DO: Need a legend still

summary.plot = function(sim.big, yvars, file_dir) {
  colors = colorRampPalette(c('red','pink','white','skyblue','darkblue'))(201)
  gamma.factor = 4 #based on the fact that gamma seems to vary btw -25 and +25, 
                   #multiplying by 4 then adding 100 gives a range of values from 0 to 200
                   #for indexing colors; may need to be adjusted if gamma exceeds this range.
  
  pdf(paste(file_dir,'/summary_plots_alpha_',Alpha,'_beta_',Beta,'.pdf',sep=''),height=6,width=8)
  par(mar = c(4,4,4,1))
  for (i in yvars) {
    y = sim.big[,which(names(sim.big)==i)]
    if (i %in% c('BK.env','BK.reg')) {
      col.index = round((y - min(y,na.rm=T))/(max(y,na.rm=T) - min(y,na.rm=T)),2)*200
    } else {
      col.index = round(y,2)*100 + 100
    }
    plot(sim.big$w.K, sim.big$sigma.reg, pch = sim.big$symbol, xlab = "<--- Environmental Filtering",
         ylab="<--- Niche Conservatism",col=colors[col.index], 
         main = paste("alpha = ",Alpha,"; beta = ",Beta,"; color = ",i,sep=''), cex=2)
    mtext("red - , blue +",3,line=0.5)
  }
  dev.off()
}

