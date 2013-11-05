#This code plots the mean +/- 2 SD of each of 4 simulation metrics:
# --global richness
# --latitude-richness correlation
# --gamma
# --MRD-richness correlation
#thru time over the course of simulations for region of origin x speciation rate

require(abind)

if(Allen==1) {
  setwd('c:/documents and settings/hurlbert/species-energy-simulation')
  sim_dir = "C:/SENCoutput/longtime_reps"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots"
} else {
  setwd('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation')
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever all of your zipped output files are
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever you want to store the results of these analyses
}

trop.sims.06 = c(4065:4072,4074)
trop.sims.07 = 4085:4094
trop.sims.08 = 4095:4104
temp.sims.06 = 4075:4084
temp.sims.07 = 4115:4124
temp.sims.08 = 4125:4134

sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv",header=T);

metric.abind = function(sims) {
  metrics = matrix(NA, nrow = 100, ncol = 28)
  for (i in sims) {
    temp = read.csv(paste(sim_dir,"/SENC_Stats_sim",i,"_mult_times.csv",sep=""),header=T)
    temp$r.lat.rich = -temp$r.env.rich
    # There is no output for timesteps in which no correlations could be calculated
    # so we add the relevant number of rows of data with NA's in that case
    if (min(temp$time) > 1000) {
      temp.top = data.frame(matrix(NA, nrow = (min(temp$time)/1000)-1, ncol = 28))
      temp.top[,3] = (1:(min(temp$time)/1000-1))*1000
      names(temp.top) = names(temp)
      temp = rbind(temp.top, temp)
    }
    metrics = abind(metrics, temp, along = 3)
  }
  return(metrics[,,-1]) #don't include the first slice of NAs
}

trop06.metrics = metric.abind(trop.sims.06)
trop07.metrics = metric.abind(trop.sims.07)
trop08.metrics = metric.abind(trop.sims.08)
temp06.metrics = metric.abind(temp.sims.06)
temp07.metrics = metric.abind(temp.sims.07)
temp08.metrics = metric.abind(temp.sims.08)

trop06.metrics.mean = data.frame(apply(trop06.metrics, 1:2, mean))
trop06.metrics.sd = data.frame(apply(trop06.metrics, 1:2, function(x) var(x)^.5))
trop07.metrics.mean = data.frame(apply(trop07.metrics, 1:2, mean))
trop07.metrics.sd = data.frame(apply(trop07.metrics, 1:2, function(x) var(x)^.5))
trop08.metrics.mean = data.frame(apply(trop08.metrics, 1:2, mean))
trop08.metrics.sd = data.frame(apply(trop08.metrics, 1:2, function(x) var(x)^.5))
temp06.metrics.mean = data.frame(apply(temp06.metrics, 1:2, mean))
temp06.metrics.sd = data.frame(apply(temp06.metrics, 1:2, function(x) var(x)^.5))
temp07.metrics.mean = data.frame(apply(temp07.metrics, 1:2, mean))
temp07.metrics.sd = data.frame(apply(temp07.metrics, 1:2, function(x) var(x)^.5))
temp08.metrics.mean = data.frame(apply(temp08.metrics, 1:2, mean))
temp08.metrics.sd = data.frame(apply(temp08.metrics, 1:2, function(x) var(x)^.5))





# Plot 4 metrics over the course of the simulation: global richness, the latitude-richness correlation, 
# gamma, and the MRD-richness correlation. Means +/- 2 SD are shown.
pdf(paste(analysis_dir,'/metrics_thru_time_varying_specn_',Sys.Date(),'.pdf',sep=""), height = 6, width = 8)
par(mfrow = c(2, 2), mar = c(3, 6, 1, 1), oma = c(4, 0, 3, 0), cex.lab = 2, las = 1, cex.axis = 1.3, mgp = c(4,1,0))
y.labels = c('Global richness', expression(italic(r)[latitude-richness]), expression(gamma), expression(italic(r)[MRD-richness]))
metric.names = c('global.richness','r.lat.rich', 'gamma.stat','r.MRD.rich')
for (j in 1:4) {
  curr.metric = metric.names[j]
  plot(trop06.metrics.mean$time/1000, trop06.metrics.mean[, curr.metric], xlim = c(0, max(trop06.metrics.mean$time/1000)), 
       ylim = range(c(trop06.metrics[, curr.metric, ], trop08.metrics[, curr.metric, ]), na.rm= T), type = "n",
       ylab = y.labels[j], xlab = "")
  polygon(c(trop06.metrics.mean$time/1000, rev(trop06.metrics.mean$time/1000)), 
          c(trop06.metrics.mean[, curr.metric] - 2*trop06.metrics.sd[, curr.metric], 
            rev(trop06.metrics.mean[, curr.metric] + 2*trop06.metrics.sd[, curr.metric])), 
          col = rgb(.8, 0, 0, .3), border = NA)
  polygon(c(trop07.metrics.mean$time/1000, rev(trop07.metrics.mean$time/1000)), 
          c(trop07.metrics.mean[, curr.metric] - 2*trop07.metrics.sd[, curr.metric], 
            rev(trop07.metrics.mean[, curr.metric] + 2*trop07.metrics.sd[, curr.metric])), 
          col = rgb(0, 0.8, 0, .3), border = NA)  
  polygon(c(trop08.metrics.mean$time/1000, rev(trop08.metrics.mean$time/1000)), 
          c(trop08.metrics.mean[, curr.metric] - 2*trop08.metrics.sd[, curr.metric], 
            rev(trop08.metrics.mean[, curr.metric] + 2*trop08.metrics.sd[, curr.metric])), 
          col = rgb(0, 0, 0.8, .3), border = NA)
  points(trop06.metrics.mean$time/1000, trop06.metrics.mean[, curr.metric], type = 'l', col = 'red', lwd = 3)
  points(trop07.metrics.mean$time/1000, trop07.metrics.mean[, curr.metric], type = 'l', col = 'green', lwd = 3)
  points(trop08.metrics.mean$time/1000, trop08.metrics.mean[, curr.metric], type = 'l', col = 'blue', lwd = 3)
  if(curr.metric == 'gamma.stat') { abline(h = 0, lty = 'dashed')}
} 
mtext("Time (x1000)", 1, outer=T, cex = 1.75, line = 2)
mtext("Tropical origin", 3, outer = T, cex = 1.5, line = 1)

# Temperate origin
for (j in 1:4) {
  curr.metric = metric.names[j]
  plot(temp06.metrics.mean$time/1000, temp06.metrics.mean[, curr.metric], xlim = c(0, max(temp06.metrics.mean$time/1000)), 
       ylim = range(c(temp06.metrics[, curr.metric, ], temp08.metrics[, curr.metric, ]), na.rm= T), type = "n",
       ylab = y.labels[j], xlab = "")
  polygon(c(temp06.metrics.mean$time/1000, rev(temp06.metrics.mean$time/1000)), 
          c(temp06.metrics.mean[, curr.metric] - 2*temp06.metrics.sd[, curr.metric], 
            rev(temp06.metrics.mean[, curr.metric] + 2*temp06.metrics.sd[, curr.metric])), 
          col = rgb(.8, 0, 0, .3), border = NA)
  polygon(c(temp07.metrics.mean$time/1000, rev(temp07.metrics.mean$time/1000)), 
          c(temp07.metrics.mean[, curr.metric] - 2*temp07.metrics.sd[, curr.metric], 
            rev(temp07.metrics.mean[, curr.metric] + 2*temp07.metrics.sd[, curr.metric])), 
          col = rgb(0, 0.8, 0, .3), border = NA)  
  polygon(c(temp08.metrics.mean$time/1000, rev(temp08.metrics.mean$time/1000)), 
          c(temp08.metrics.mean[, curr.metric] - 2*temp08.metrics.sd[, curr.metric], 
            rev(temp08.metrics.mean[, curr.metric] + 2*temp08.metrics.sd[, curr.metric])), 
          col = rgb(0, 0, 0.8, .3), border = NA)
  points(temp06.metrics.mean$time/1000, temp06.metrics.mean[, curr.metric], type = 'l', col = 'red', lwd = 3)
  points(temp07.metrics.mean$time/1000, temp07.metrics.mean[, curr.metric], type = 'l', col = 'green', lwd = 3)
  points(temp08.metrics.mean$time/1000, temp08.metrics.mean[, curr.metric], type = 'l', col = 'blue', lwd = 3)
  if(curr.metric == 'gamma.stat') { abline(h = 0, lty = 'dashed')}
} 
mtext("Time (x1000)", 1, outer=T, cex = 1.75, line = 2)
mtext("Temperate origin", 3, outer = T, cex = 1.5, line = 1)

dev.off()
