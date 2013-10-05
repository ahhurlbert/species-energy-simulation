## Showing how various simulation metrics vary through time and as a function
# of sigma_E (low sigma = strong niche conservatism, high sigma = weak niche conservatism)
require(abind)
require(apTreeshape)


Allen = 0;

#New parameter for taking into account which of us is running this code
if(Allen==1) {
  setwd('c:/documents and settings/hurlbert/species-energy-simulation')
  sim_dir = "C:/SENCoutput"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots"
} else {
  setwd('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation')
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever all of your zipped output files are
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever you want to store the results of these analyses
}

#Energy gradient sims, 10 sims each

#sigma_E = 1
Ktrop.sims.sig1 = 4065:4074
Ktemp.sims.sig1 = 4075:4084

#sigma_E = 3
Ktrop.sims.sig3 = 2635:2644
Ktemp.sims.sig3 = 2795:2804

#sigma_E = 9
Ktrop.sims.sig9 = 2675:2684
Ktemp.sims.sig9 = 2835:2844

sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv",header=T);

metric.abind.new = function(sims, scenario = "K", min.div.regions = 4, min.richness = 30) {
  num.cols = 41
  metrics = matrix(NA, nrow = 100, ncol = num.cols)
  for (i in sims) {
    if (i < 4065) {
      temp = read.csv(paste(sim_dir,"/NEW_Stats_sim",i,"_time_seq_root_only.csv",sep=""),header=T)
    } else if (i >= 4065) {
      temp = read.csv(paste(sim_dir,"/lbeta_Stats_sim",i,"_time_seq_root_only.csv",sep=""),header=T)
    }
    temp$r.lat.rich = -temp$r.env.rich
    temp$scaled.MRD.range = temp$MRD.range/temp$max.RD
    temp$scaled.MRD.rich.slope = temp$MRD.rich.slope/temp$max.RD
    # There is no output for timesteps in which no correlations could be calculated
    # so we add the relevant number of rows of data with NA's in that case
    if (nrow(temp) < 100) {
      temp.top = data.frame(matrix(NA, nrow = 100 - nrow(temp), ncol = num.cols))
      names(temp.top) = names(temp)
      temp = rbind(temp.top, temp)
    }
    # Only include metrics in analysis if simulation meets minimum region and richness criteria
    if (min(temp$n.div.regions, na.rm = T) < min.div.regions) {
      temp[which(temp$n.div.regions < min.div.regions),] = NA
    }
    if (min(temp$global.richness, na.rm = T) < min.richness) {
      temp[which(temp$global.richness < min.richness),] = NA
    }
    metrics = abind(metrics, temp, along = 3)
    #print(c(i,range(temp$global.richness,na.rm=T),range(temp$n.div.regions,na.rm=T)))
  }
  return(metrics[,,-1]) #don't include the first slice of NAs
}

min.num.regions = 5
min.num.div.regions = 5
min.global.richness = 30

Ktemp.sig1.metrics = metric.abind.new(Ktemp.sims.sig1, scenario = "K", min.div.regions = min.num.div.regions, min.richness = min.global.richness)
Ktrop.sig1.metrics = metric.abind.new(Ktrop.sims.sig1, scenario = "K", min.div.regions = min.num.div.regions, min.richness = min.global.richness)
Ktemp.sig3.metrics = metric.abind.new(Ktemp.sims.sig3, scenario = "K", min.div.regions = min.num.div.regions, min.richness = min.global.richness)
Ktrop.sig3.metrics = metric.abind.new(Ktrop.sims.sig3, scenario = "T", min.div.regions = min.num.div.regions, min.richness = min.global.richness)
Ktemp.sig9.metrics = metric.abind.new(Ktemp.sims.sig9, scenario = "K", min.div.regions = min.num.div.regions, min.richness = min.global.richness)
Ktrop.sig9.metrics = metric.abind.new(Ktrop.sims.sig9, scenario = "T", min.div.regions = min.num.div.regions, min.richness = min.global.richness)

#Function for calculating mean or SD for simulations with a minimum number of non-NA values at a given time step
calc.meanSD = function(x, stat = 'mean', min.num.nonNA = 10) {
  if (stat == 'mean') {
    if(sum(!is.na(x)) >= min.num.nonNA) {
      mean(x, na.rm = T)
    } else { NA }
  } else if (stat == 'sd') {
    if(sum(!is.na(x)) >= min.num.nonNA) {
      var(x, na.rm = T)^0.5
    } else { NA }
  }
}

min.num.datapts = 10

Ktemp.sig1.mean = data.frame(apply(Ktemp.sig1.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
Ktemp.sig1.sd = data.frame(apply(Ktemp.sig1.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))

Ktrop.sig1.mean = data.frame(apply(Ktrop.sig1.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
Ktrop.sig1.sd = data.frame(apply(Ktrop.sig1.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))

Ktemp.sig3.mean = data.frame(apply(Ktemp.sig3.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
Ktemp.sig3.sd = data.frame(apply(Ktemp.sig3.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))

Ktrop.sig3.mean = data.frame(apply(Ktrop.sig3.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
Ktrop.sig3.sd = data.frame(apply(Ktrop.sig3.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))

Ktemp.sig9.mean = data.frame(apply(Ktemp.sig9.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
Ktemp.sig9.sd = data.frame(apply(Ktemp.sig9.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))

Ktrop.sig9.mean = data.frame(apply(Ktrop.sig9.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
Ktrop.sig9.sd = data.frame(apply(Ktrop.sig9.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))

metric.names = c('global.richness',
                 'r.lat.rich', 
                 'gamma.stat',
                 'r.env.PSV', 
                 'r.env.MRD', 
                 'r.MRD.rich',
                 'r.PSV.rich',
                 'MRD.rich.slope',
                 'MRD.env.slope',
                 'PSV.rich.slope',
                 'PSV.env.slope',
                 'MRD.range',
                 'PSV.range',
                 'MRD.mean',
                 'PSV.mean',
                 'tree.beta',
                 'scaled.MRD.range',
                 'scaled.MRD.rich.slope')

metric.labels = c('Global richness', 
                  expression(italic(r)[latitude-richness]), 
                  expression(gamma), 
                  expression(italic(r)[env-PSV]),
                  expression(italic(r)[env-MRD]), 
                  expression(italic(r)[MRD-richness]),
                  expression(italic(r)[PSV-richness]),
                  'MRD-Richness slope',
                  'MRD-Environment slope',
                  'PSV-Richness slope',
                  'PSV-Environment slope',
                  'MRD range',
                  'PSV range',
                  'Mean MRD',
                  'Mean PSV',
                  expression(beta),
                  'scaled MRD range',
                  'MRD-richness slope')


# Plotting metrics over the course of the simulation: EXPLORATORY PLOTS
pdf(paste(analysis_dir, '/metrics_thru_time_sigma_dependence', Sys.Date(), '.pdf', sep = ""), height = 12, width = 10)
par(mfrow = c(3, 2), mar = c(5, 6, 1, 1), oma = c(5, 0, 0, 0), cex.lab = 2, las = 1, cex.axis = 1.3, mgp = c(4,1,0))

# Specify variables to plot here, and width of error bars
names4plotting = c('global.richness','r.lat.rich', 'scaled.MRD.rich.slope', 'tree.beta', 'PSV.rich.slope', 'gamma.stat')

error = 2 # error bars in SD units (+/-)
for (j in 1:6) {
  curr.metric = names4plotting[j]

  #Sigma_E = 1
  plot(Ktrop.sig1.mean$time/1000, Ktrop.sig1.mean[, curr.metric], xlim = c(0, max(Ktrop.sig1.mean$time, na.rm=T)/1000), 
       ylim = range(c( Ktemp.sig1.metrics[, curr.metric, ], 
                       Ktemp.sig3.metrics[, curr.metric, ],
                       Ktemp.sig9.metrics[, curr.metric, ]),na.rm= T), 
       type = "n",ylab = metric.labels[metric.names == curr.metric], xlab = "")
  #polygon(c(Ktrop.sig1.mean$time/1000, rev(Ktrop.sig1.mean$time/1000)), 
  #        c(Ktrop.sig1.mean[, curr.metric] - error*Ktrop.sig1.sd[, curr.metric], 
  #          rev(Ktrop.sig1.mean[, curr.metric] + error*Ktrop.sig1.sd[, curr.metric])), 
  #        col = rgb(.8, 0.8, 0.4, .3), border = NA)
  polygon(c(Ktemp.sig1.mean$time/1000, rev(Ktemp.sig1.mean$time/1000)), 
          c(Ktemp.sig1.mean[, curr.metric] - error*Ktemp.sig1.sd[, curr.metric], 
            rev(Ktemp.sig1.mean[, curr.metric] + error*Ktemp.sig1.sd[, curr.metric])), 
          col = rgb(0.8, 0.8, 0, .3), border = NA)
  #points(Ktrop.sig1.mean$time/1000, Ktrop.sig1.mean[, curr.metric], type = 'l', col = rgb(1, 1, .5), lwd = 3)
  points(Ktemp.sig1.mean$time/1000, Ktemp.sig1.mean[, curr.metric], type = 'l', col = rgb(1, 1, .5), lwd = 3)
  
  par(new = T)
  
  #sigma_E = 3
  plot(Ktrop.sig3.mean$time/1000, Ktrop.sig3.mean[, curr.metric], xlim = c(0, max(Ktrop.sig3.mean$time, na.rm=T)/1000), 
       ylim = range(c( Ktemp.sig1.metrics[, curr.metric, ], 
                       Ktemp.sig3.metrics[, curr.metric, ],
                     Ktemp.sig9.metrics[, curr.metric, ]),na.rm= T), 
       type = "n",ylab = metric.labels[metric.names == curr.metric], xlab = "", xaxt = "n", yaxt = "n")
  #polygon(c(Ktrop.sig3.mean$time/1000, rev(Ktrop.sig3.mean$time/1000)), 
  #        c(Ktrop.sig3.mean[, curr.metric] - error*Ktrop.sig3.sd[, curr.metric], 
  #          rev(Ktrop.sig3.mean[, curr.metric] + error*Ktrop.sig3.sd[, curr.metric])), 
  #        col = rgb(0.8, 0.4, 0.8, .3), border = NA)
  polygon(c(Ktemp.sig3.mean$time/1000, rev(Ktemp.sig3.mean$time/1000)), 
          c(Ktemp.sig3.mean[, curr.metric] - error*Ktemp.sig3.sd[, curr.metric], 
            rev(Ktemp.sig3.mean[, curr.metric] + error*Ktemp.sig3.sd[, curr.metric])), 
          col = rgb(0.8, 0, 0.8, .3), border = NA)
  #points(Ktrop.sig3.mean$time/1000, Ktrop.sig3.mean[, curr.metric], type = 'l', col = rgb(1, .5, 1), lwd = 3)
  points(Ktemp.sig3.mean$time/1000, Ktemp.sig3.mean[, curr.metric], type = 'l', col = rgb(1, .5, 1), lwd = 3)
  
  par(new = T)
  
  #sigma_E = 9
  plot(Ktrop.sig9.mean$time/1000, Ktrop.sig9.mean[, curr.metric], xlim = c(0, max(Ktrop.sig9.mean$time, na.rm=T)/1000), 
       ylim = range(c( Ktemp.sig1.metrics[, curr.metric, ], 
                       Ktemp.sig3.metrics[, curr.metric, ],
                       Ktemp.sig9.metrics[, curr.metric, ]),na.rm= T), 
       type = "n",ylab = metric.labels[metric.names == curr.metric], xlab = "", xaxt = "n", yaxt = "n")
  #polygon(c(Ktrop.sig9.mean$time/1000, rev(Ktrop.sig9.mean$time/1000)), 
  #        c(Ktrop.sig9.mean[, curr.metric] - error*Ktrop.sig9.sd[, curr.metric], 
  #          rev(Ktrop.sig9.mean[, curr.metric] + error*Ktrop.sig9.sd[, curr.metric])), 
  #        col = rgb(0.4, 0.8, 0.8, .3), border = NA)
  polygon(c(Ktemp.sig9.mean$time/1000, rev(Ktemp.sig9.mean$time/1000)), 
          c(Ktemp.sig9.mean[, curr.metric] - error*Ktemp.sig9.sd[, curr.metric], 
            rev(Ktemp.sig9.mean[, curr.metric] + error*Ktemp.sig9.sd[, curr.metric])), 
          col = rgb(0, 0.8, 0.8, .3), border = NA)
  #points(Ktrop.sig9.mean$time/1000, Ktrop.sig9.mean[, curr.metric], type = 'l', col = rgb(.5, 1, 1), lwd = 3)
  points(Ktemp.sig9.mean$time/1000, Ktemp.sig9.mean[, curr.metric], type = 'l', col = rgb(.5, 1, 1), lwd = 3)
  
  if (curr.metric == 'gamma.stat') { abline(h = 0, lty = 'dashed')}
  if (curr.metric == 'r.lat.rich') {
    legend('topright', legend = paste('sigma =', c(1,3,9)), col = c(rgb(1,1,.5), rgb(1,.5,1), rgb(.5,1,1)), lwd = 3, cex = 2)
  }
} 
mtext("Time (x1000)", 1, outer=T, cex = 1.75, line = 1.5) 
dev.off()



