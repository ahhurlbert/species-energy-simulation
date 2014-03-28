# Function for calculating origination and extinction rates and appending them
# as new columns to old stats output files (in which stats are calculated for
# the root clade over different timeslices)

#source('code/supplemental_analysis_functions.r')

add.orig.ext.rates = function(sim, timeslices) {
  if (file.access(paste('raw_sim_output/sim', sim, '_out/SENC_all.pops_sim', sim, '.csv', sep = ''), mode = 0) == 0) { 
    all.populations = read.csv(paste('raw_sim_output/sim', sim, '_out/SENC_all.pops_sim', sim, '.csv', sep = ''), header=T)
  } else if( file.access(paste(output_dir,'/senc.out.', sim, '.zip', sep = ''), mode = 0) == 0 ) {
    unzipped.filenames = unzip(paste('raw_sim_output/senc.out.',sim,'.zip',sep=''), 
                               exdir = paste('raw_sim_output/sim', sim, '_out', sep = ''))
    all.populations = read.csv(unzipped.filenames[1], header=T)
  }
  all.populations = read.csv(paste("raw_sim_output/sim", sim, "_out/SENC_all.pops_sim", sim, ".csv", sep = ""))
  stats = read.csv(paste("analysis_output/Stats_sim", sim, "_rootclade_only_", length(timeslices), "_times.csv", sep = ""))
  
  glob.pops = unique(all.populations[, c('spp.name', 'time.of.sp.origin', 'time.of.sp.extinction')])
  
  for (t in timeslices) {
    t0 = max(c(0, timeslices[timeslices < t]))
    if (num.of.time.slices == 1) {
      origins = length(glob.pops$spp.name[glob.pops$time.of.sp.extinction > t])
      extincts = length(glob.pops$spp.name[glob.pops$time.of.sp.extinction < t])
    } else {
      origin.pops = subset(glob.pops, time.of.sp.origin >= t0 &
                             time.of.sp.origin < t &
                             time.of.sp.extinction > t)
      origins = nrow(origin.pops)
      extinct.pops = subset(glob.pops, time.of.sp.origin <= t0 &
                              time.of.sp.extinction >= t0 &
                              time.of.sp.extinction < t)
      extincts = nrow(extinct.pops)
    }
    time.interval = t - t0
    richness.t0 = nrow( subset(glob.pops, time.of.sp.origin <= t0 & time.of.sp.extinction > t0) )
    
    output = rbind(output, cbind(glob.origins = origins, 
                                 glob.extinctions = extincts,
                                 glob.orig.rate = origins/time.interval/richness.t0,
                                 glob.ext.rate = extincts/time.interval/richness.t0) )
    
  }
  output = data.frame(output)
  stats.out = cbind(stats, output)
}