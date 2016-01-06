# Plotting extinction and origination rates thru time for
# 4 diversification scenarios presented in Hurlbert & Stegen 2014
# Frontiers in Genetics

setwd('z:/Manuscripts/FrontiersTropicalDiversity/analysis_output')

metric.abind.new = function(sims, num.cols = 45, num.rows = 100) {
  require(abind)
  metrics = matrix(NA, nrow = num.rows, ncol = num.cols)
  for (i in sims) {
    statsfile = list.files(pattern = paste("Stats_sim", i, "_rootclade", sep = ""))
    temp = read.csv(statsfile, header = T)
    
    # There is no output for timesteps in which no correlations could be calculated
    # so we add the relevant number of rows of data with NA's in that case
    if (nrow(temp) < num.rows) {
      temp.top = data.frame(matrix(NA, nrow = num.rows - nrow(temp), ncol = num.cols))
      names(temp.top) = names(temp)
      temp = rbind(temp.top, temp)
    }
    metrics = abind(metrics, temp, along = 3)
  }
  return(metrics[,,-1]) #don't include the first slice of NAs
}

# Function to plot origination and extinction rates over time for a set of scenarios
orig.ext.thru.time = function(simresults, max.time) {
  spec.color = rgb(255, 192, 0, maxColorValue = 255)
  ext.color = rgb(0, 176, 240, maxColorValue = 255)
  max.rate = max(c(simresults[, 'glob.orig.rate',], simresults[, 'glob.ext.rate',]), na.rm = T)
  min.rate = min(c(simresults[, 'glob.orig.rate',], simresults[, 'glob.ext.rate',]), na.rm = T)
  
  if (min.rate == 0) { min.rate = 1e-6 }
  plot(simresults[, 'time', 1], log10(simresults[, 'glob.orig.rate', 1]), type = 'l',
       xlab = '', ylab = '', ylim = log10(c(min.rate, max.rate)), xlim = c(0, max.time),
       col = spec.color, lwd = 2, cex.axis = 1.5)
  for (i in 2:10) {
    points(simresults[, 'time', i], log10(simresults[, 'glob.orig.rate', i]), type = 'l',
           col = spec.color, lwd = 2)
  }
  for (i in 1:10) {
    points(simresults[, 'time', i], log10(simresults[, 'glob.ext.rate', i]), type = 'l', 
           col = ext.color, lwd = 2)
  }
}

#Niche conservatism sims
ncsims = metric.abind.new(3465:3474, num.cols = 42)
#Energy gradient sims
egsims = metric.abind.new(4065:4074, num.cols = 42)
#Speciation gradient sims
sgsims = metric.abind.new(5525:5534, num.cols = 42)
#Disturbance gradient sims
dgsims = metric.abind.new(5625:5634, num.cols = 42)

pdf('z:/git/species-energy-simulation/analysis_output/4scenarios_rates_thru_time.pdf',
    height = 4, width = 12)
par(mfrow = c(1, 4), oma = c(4, 4, 0, 0), mar = c(3, 3, 1, 1), las = 1)
orig.ext.thru.time(egsims, max.time = 25000)
orig.ext.thru.time(sgsims, max.time = 25000)
orig.ext.thru.time(dgsims, max.time = 25000)
orig.ext.thru.time(ncsims, max.time = 200)
mtext("Time", 1, outer = T, cex = 2, line = 2)
mtext("log10 Rate", 2, outer = T, cex = 2, line = 1, las = 0)
dev.off()

