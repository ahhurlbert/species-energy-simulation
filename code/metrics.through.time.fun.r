## code for relating various metrics of the root clade to time and making plots

# Function for abinding simulation output (the "stats" files where rows are points in time)
# across all of the replicates with identical parameters
# --current stats output is 42 cols, plus 3 more cols which are created within the function
# --older stats output may have fewer cols

metric.abind.new = function(sims, min.div.regions = 4, min.richness = 30, num.cols = 45) {
  require(abind)
  metrics = matrix(NA, nrow = 100, ncol = num.cols)
  for (i in sims) {
    statsfile = list.files(path = "analysis_output", pattern = paste("Stats_sim", i, "_rootclade", sep = ""))
    temp = read.csv(paste("analysis_output/", statsfile, sep = ""), header = T)
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
  }
  return(metrics[,,-1]) #don't include the first slice of NAs
}


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


plot.metrics.thru.time = function(trop.sims, 
                                  temp.sims, 
                                  sim.matrix,
                                  pdf.out,
                                  min.div.regions = 4,
                                  min.richness = 30,
                                  min.num.datapts = 10,
                                  num.cols = 45)
{
  #require(apTreeshape)
  #require(ape)

  temp.metrics = metric.abind.new(temp.sims, min.div.regions = min.div.regions, min.richness = min.richness, num.cols = num.cols)
  trop.metrics = metric.abind.new(trop.sims, min.div.regions = min.div.regions, min.richness = min.richness, num.cols = num.cols)
  
  temp.metrics.mean = data.frame(apply(temp.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
  temp.metrics.sd = data.frame(apply(temp.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))
  
  trop.metrics.mean = data.frame(apply(trop.metrics, 1:2, function(x) calc.meanSD(x, stat = 'mean', min.num.nonNA = min.num.datapts)))
  trop.metrics.sd = data.frame(apply(trop.metrics, 1:2, function(x) calc.meanSD(x, stat = 'sd', min.num.nonNA = min.num.datapts)))
  
  metric.names = c('global.richness',
                   'r.lat.rich', 
                   'gamma.stat',
                   'scaled.MRD.rich.slope',
                   'tree.beta')
  
  metric.labels = c('Global richness', 
                    expression(italic(r)[latitude-richness]), 
                    expression(gamma), 
                    'MRD-Richness slope',
                    expression(beta))

  pdf(paste('analysis_output/', pdf.out, sep = ''), width = 11, height = 8)
  par(mfrow = c(2, 3), mar = c(5, 6, 1, 1), oma = c(5, 0, 5, 0), cex.lab = 2, las = 1, cex.axis = 1.3, mgp = c(4,1,0))

  sim.params = sim.matrix[sim.matrix$sim.id == trop.sims[1], ]
  if (sim.params$disturb_frequency == 0) {disturb = 'no'} else {disturb = 'yes'}
  if (sim.params$specn.gradient == 'off') {spec.grad = 'off'} else {spec.grad = paste(sim.params$specn.factor, "x", sep = "")}
  if (sim.params$carry.cap == 'off') { 
    x.offset.tmp = min(temp.metrics.mean$time, na.rm = T)
    x.offset.trp = min(trop.metrics.mean$time, na.rm = T)
    x.scale = 1
  } else {
    x.offset.tmp = 0
    x.offset.trp = 0
    x.scale = 1000
  }
  
  error = 2 # error bars in SD units (+/-)
  for (curr.metric in metric.names) {
    plot(trop.metrics.mean$time/x.scale - x.offset.tmp, trop.metrics.mean[, curr.metric], 
         xlim = c(0, max(trop.metrics.mean$time, na.rm=T)/x.scale - x.offset.tmp), 
         ylim = range(c(trop.metrics[, curr.metric, ], temp.metrics[, curr.metric, ]), na.rm= T), type = "n",
         ylab = metric.labels[metric.names == curr.metric], xlab = "")
    polygon(c(trop.metrics.mean$time/x.scale - x.offset.trp, rev(trop.metrics.mean$time/x.scale - x.offset.trp)), 
            c(trop.metrics.mean[, curr.metric] - error*trop.metrics.sd[, curr.metric], 
              rev(trop.metrics.mean[, curr.metric] + error*trop.metrics.sd[, curr.metric])), 
            col = rgb(.8, 0, 0, .3), border = NA)
    polygon(c(temp.metrics.mean$time/x.scale - x.offset.tmp, rev(temp.metrics.mean$time/x.scale - x.offset.tmp)), 
            c(temp.metrics.mean[, curr.metric] - error*temp.metrics.sd[, curr.metric], 
              rev(temp.metrics.mean[, curr.metric] + error*temp.metrics.sd[, curr.metric])), 
            col = rgb(0, 0, .8, .3), border = NA)
    points(trop.metrics.mean$time/x.scale - x.offset.trp, trop.metrics.mean[, curr.metric], type = 'l', col = 'red', lwd = 3)
    points(temp.metrics.mean$time/x.scale - x.offset.tmp, temp.metrics.mean[, curr.metric], type = 'l', col = 'blue', lwd = 3)
    
    if(curr.metric != 'global.richness') { abline(h = 0, lty = 'dashed')}
  }
  # extinction and origination rates panel
  rate.range = log10(range(c(temp.metrics.mean$glob.orig.rate[temp.metrics.mean$glob.orig.rate > 0 & 
                                                                is.na(temp.metrics.mean$glob.orig.rate) == F], 
                             temp.metrics.mean$glob.ext.rate[temp.metrics.mean$glob.ext.rate > 0 & 
                                                               is.na(temp.metrics.mean$glob.orig.rate) == F])))
  plot(temp.metrics.mean$time/x.scale - x.offset.tmp, log10(temp.metrics.mean$glob.orig.rate), type = 'l', col = 'blue', lwd =3, 
       xlab = "", ylab = "Rate", ylim = rate.range)
  points(temp.metrics.mean$time/x.scale - x.offset.tmp, log10(temp.metrics.mean$glob.ext.rate), type = 'l', col = 'blue', lwd =3, lty = 'dashed')
  points(trop.metrics.mean$time/x.scale - x.offset.trp, log10(trop.metrics.mean$glob.orig.rate), type = 'l', col = 'red', lwd =3)
  points(trop.metrics.mean$time/x.scale - x.offset.trp, log10(trop.metrics.mean$glob.ext.rate), type = 'l', col = 'red', lwd =3, lty = 'dashed')
  legend("topright", c('origination rate', 'extinction rate'), lty = c('solid', 'dashed'))
  
  mtext("Time", 1, outer=T, cex = 1.75, line = 1.5) 
  mtext(paste("Sims", min(c(temp.sims, trop.sims)), "-", max(c(temp.sims, trop.sims)), "; Energetic constraint", sim.params$carry.cap[1], "; K gradient", sim.params$energy.gradient[1], "; w =",
        sim.params$w[1], ";\ngamma =", sim.params$gamma[1], "; sigma =", sim.params$sigma_E[1], "; disturbance =", 
        disturb, "; speciation gradient", spec.grad,
        ";\nalpha =", sim.params$alpha[1], "; disp =", sim.params$beta[1]), 3, outer=T, line = .5)
  dev.off()
}