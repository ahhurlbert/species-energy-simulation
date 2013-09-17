## code for relating various metrics of the root clade to time and making plots
require(abind)



Allen = 0;

#New parameter for taking into account which of us is running this code
if(Allen==1) {
  setwd('c:/documents and settings/hurlbert/species-energy-simulation')
  sim_dir = "C:/SENCoutput/longtime_reps"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots"
} else {
  setwd('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation')
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever all of your zipped output files are
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever you want to store the results of these analyses
}

#Energy gradient sims
trop.sims = 4065:4074
temp.sims = 4075:4084

#No energy gradient sims
Ttrop.sims = 3465:3474
Ttemp.sims = 3565:3574

sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv",header=T);

metric.abind = function(sims, scenario = "K") {
  metrics = matrix(NA, nrow = 100, ncol = 28)
  for (i in sims) {
    if (scenario == "K") {
      temp = read.csv(paste(sim_dir,"/SENC_Stats_sim",i,"_mult_times.csv",sep=""),header=T)
    } else if (scenario == "T") {
      temp = read.csv(paste(sim_dir,"/SENC_Stats_sim",i,"_mult_times_root_only.csv",sep=""),header=T)
    }
    temp$r.lat.rich = -temp$r.env.rich
    # There is no output for timesteps in which no correlations could be calculated
    # so we add the relevant number of rows of data with NA's in that case
    if (nrow(temp) < 100) {
      temp.top = data.frame(matrix(NA, nrow = 100 - nrow(temp), ncol = 28))
      names(temp.top) = names(temp)
      temp = rbind(temp.top, temp)
    }
    metrics = abind(metrics, temp, along = 3)
  }
  return(metrics[,,-1]) #don't include the first slice of NAs
}

temp.metrics = metric.abind(temp.sims, scenario = "K")
trop.metrics = metric.abind(trop.sims, scenario = "K")
Ttemp.metrics = metric.abind(Ttemp.sims, scenario = "T")
Ttrop.metrics = metric.abind(Ttrop.sims, scenario = "T")


temp.metrics.mean = data.frame(apply(temp.metrics, 1:2, mean))
temp.metrics.sd = data.frame(apply(temp.metrics, 1:2, function(x) var(x)^.5))

trop.metrics.mean = data.frame(apply(trop.metrics, 1:2, mean))
trop.metrics.sd = data.frame(apply(trop.metrics, 1:2, function(x) var(x)^.5))

Ttemp.metrics.mean = data.frame(apply(Ttemp.metrics, 1:2, mean))
Ttemp.metrics.sd = data.frame(apply(Ttemp.metrics, 1:2, function(x) var(x)^.5))

Ttrop.metrics.mean = data.frame(apply(Ttrop.metrics, 1:2, mean))
Ttrop.metrics.sd = data.frame(apply(Ttrop.metrics, 1:2, function(x) var(x)^.5))

# Plot 4 metrics over the course of the simulation: global richness, the latitude-richness correlation, 
# gamma, and the MRD-richness correlation. Means +/- 2 SD are shown.
pdf(paste(analysis_dir,'/metrics_thru_time_inc_Tscenario',Sys.Date(),'.pdf',sep=""), height = 6, width = 8)
par(mfrow = c(2, 2), mar = c(3, 6, 1, 1), oma = c(3, 0, 0, 0), cex.lab = 2, las = 1, cex.axis = 1.3, mgp = c(4,1,0))
metric.names = c('global.richness','r.lat.rich', 'gamma.stat','r.env.PSV', 'r.env.MRD', 'r.MRD.rich','r.PSV.rich')
metric.labels = c('Global richness', expression(italic(r)[latitude-richness]), 
                  expression(gamma), expression(italic(r)[env-PSV]),
                  expression(italic(r)[env-MRD]), expression(italic(r)[MRD-richness]),
                  expression(italic(r)[PSV-richness]))

# Specify variables to plot here, and width of error bars
names4plotting = c('global.richness','r.lat.rich', 'gamma.stat','r.MRD.rich')
error = 2 # error bars in SD units (+/-)
for (j in 1:4) {
  curr.metric = names4plotting[j]
  plot(trop.metrics.mean$time/1000, trop.metrics.mean[, curr.metric], xlim = c(0, max(trop.metrics.mean$time/1000)), 
       ylim = range(c(trop.metrics[, curr.metric, ], temp.metrics[, curr.metric, ]), na.rm= T), type = "n",
       ylab = metric.labels[metric.names == curr.metric], xlab = "")
  polygon(c(trop.metrics.mean$time/1000, rev(trop.metrics.mean$time/1000)), 
          c(trop.metrics.mean[, curr.metric] - error*trop.metrics.sd[, curr.metric], 
            rev(trop.metrics.mean[, curr.metric] + error*trop.metrics.sd[, curr.metric])), 
          col = rgb(.8, 0, 0, .3), border = NA)
  polygon(c(temp.metrics.mean$time/1000, rev(temp.metrics.mean$time/1000)), 
          c(temp.metrics.mean[, curr.metric] - error*temp.metrics.sd[, curr.metric], 
            rev(temp.metrics.mean[, curr.metric] + error*temp.metrics.sd[, curr.metric])), 
          col = rgb(0, 0, .8, .3), border = NA)
  points(trop.metrics.mean$time/1000, trop.metrics.mean[, curr.metric], type = 'l', col = 'red', lwd = 3)
  points(temp.metrics.mean$time/1000, temp.metrics.mean[, curr.metric], type = 'l', col = 'blue', lwd = 3)
  
  par(new = T)
  
  plot(Ttrop.metrics.mean$time - min(Ttrop.metrics.mean$time, na.rm = T), Ttrop.metrics.mean[, curr.metric], 
       #xlim = c(0, max(Ttrop.metrics.mean$time, na.rm = T) - min(Ttrop.metrics.mean$time, na.rm = T)), 
       ylim = range(c(Ttrop.metrics[, curr.metric, ], Ttemp.metrics[, curr.metric, ]), na.rm= T), type = "n",
       ylab = "", xlab = "", yaxt = "n", xaxt = "n")
  polygon(c(Ttrop.metrics.mean$time - min(Ttrop.metrics.mean$time, na.rm = T), rev(Ttrop.metrics.mean$time - min(Ttrop.metrics.mean$time, na.rm = T))), 
          c(Ttrop.metrics.mean[, curr.metric] - error*Ttrop.metrics.sd[, curr.metric], 
            rev(Ttrop.metrics.mean[, curr.metric] + error*Ttrop.metrics.sd[, curr.metric])), 
          col = rgb(.8, 0, 0, .3), border = NA)
  polygon(c(Ttemp.metrics.mean$time - min(Ttemp.metrics.mean$time, na.rm = T), rev(Ttemp.metrics.mean$time - min(Ttemp.metrics.mean$time, na.rm = T))), 
          c(Ttemp.metrics.mean[, curr.metric] - error*Ttemp.metrics.sd[, curr.metric], 
            rev(Ttemp.metrics.mean[, curr.metric] + error*Ttemp.metrics.sd[, curr.metric])), 
          col = rgb(0, 0, .8, .3), border = NA)
  points(Ttrop.metrics.mean$time - min(Ttrop.metrics.mean$time, na.rm = T), Ttrop.metrics.mean[, curr.metric], type = 'l', col = 'red', lwd = 3, lty = 'dashed')
  points(Ttemp.metrics.mean$time - min(Ttemp.metrics.mean$time, na.rm = T), Ttemp.metrics.mean[, curr.metric], type = 'l', col = 'blue', lwd = 3, lty = 'dashed')
  
  
  if(curr.metric == 'gamma.stat') { abline(h = 0, lty = 'dashed')}
} 
mtext("Time (x1000)", 1, outer=T, cex = 1.75, line = 1.5)                                      
dev.off()
