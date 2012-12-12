# Evaluating relationship between simulation output correlations and simulation parameters

file_dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots'

#sim.rs = compile.firstlines(file_dir , "SENC_Stats_sim")          #takes ~12 minutes on Hurlbert office machine
#write.csv(sim.rs,paste(file_dir,'/allsims_bigclade_stats_output.csv',sep=''), row.names=F)

sim.rs = read.csv(paste(file_dir,'/allsims_bigclade_stats_output.csv',sep=''), header=T)
sim.matrix = read.csv(paste(file_dir,'/sim.matrix.output_2012-12-10.csv',sep=''), header=T)

sim.big = merge(sim.matrix, sim.rs, by.x="sim.id", by.y="sim", all.x=T)
sim.big = sim.big[sim.big$sim.id!=0,]

# List of independent variables to plot (using color)
yvars = c('r.time.rich','r.env.rich','r.MRD.rich','r.PSV.rich','r.env.MRD',
          'r.env.PSV','r.ext.reg','r.rich.ext')


# Plots for K gradient simulations
# -- plots are on sigma_E versus w space, with symbol size reflecting beta/alpha
# -- and color reflecting yvars (list above)

summary.plot = function(sim.big, carry.cap, energy.gradient, yvars, file_dir) {
  colors = colorRampPalette(c('red','pink','white','skyblue','darkblue'))(201)
  gamma.factor = 4 #based on the fact that gamma seems to vary btw -25 and +25, 
                   #multiplying by 4 then adding 100 gives a range of values from 0 to 200
                   #for indexing colors; may need to be adjusted if gamma exceeds this range.
  
  # Subset the data according to the arguments
  sub.sim = subset(sim.big, carry.cap==carry.cap & energy.gradient==energy.gradient)
  
  pdf(paste(file_dir,'/K_',carry.cap,'_Grad_',energy.gradient,'_summary_plots.pdf',sep=''),height=6,width=8)
  for (i in yvars) {
    plot(jitter(sub.sim$w),jitter(sub.sim$sigma_E),pch=16, xlab = "<--- Environmental Filtering",
         ylab="<--- Niche Conservatism",col=colors[(round(sub.sim[,which(names(sub.sim)==i)],2)*100)+100], 
         main = paste("K ",carry.cap,", Gradient ",energy.gradient,", size = beta/alpha, color = ",i,sep=''), 
         cex = log10(sub.sim$beta/sub.sim$alpha)*.5+.5)
    mtext("red - , blue +",3,line=0.5)
  }
  plot(jitter(sub.sim$w),jitter(sub.sim$sigma_E),pch=16, xlab = "<--- Environmental Filtering",
       ylab="<--- Niche Conservatism",col=colors[round(sub.sim$gamma.stat*gamma.factor,0)+100], 
       main = paste("K ",carry.cap,", Gradient ",energy.gradient,", size = beta/alpha, color = gamma",sep=''), 
       cex = log10(sub.sim$beta/sub.sim$alpha)*.5+.5)
  mtext("red - , blue +",3,line=0.5)
  dev.off()
}

