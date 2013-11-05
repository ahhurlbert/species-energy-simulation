# Script for visualizing output across simulations for the entire clade
# (i.e., the first row of the sim output for each sim #)

# First, get a list of the proper files, pull out the first line of each and store together in a df
file_dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses'
filetext = "SENC_Stats_sim"

compile.firstlines = function(file_dir, filetext) {
  files = list.files(file_dir)
  stats.files = files[grep(filetext,files)]
  simstats.finaltime = c()
  for (i in stats.files) {
    stats = read.csv(paste(file_dir,'/',i,sep=''), header=T)
    if (nrow(stats) > 0) {
      simstats.finaltime = rbind(simstats.finaltime, stats[1,])
    } else {
      sim = as.numeric( strsplit( strsplit(i,"m")[[1]][2] , "\\.")[[1]][1] )
      simstats.finaltime = rbind(simstats.finaltime, c(sim, rep(NA,(ncol(simstats.finaltime)-1))))
    }
  }
  return(simstats.finaltime)
}

