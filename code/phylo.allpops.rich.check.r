# Function for checking simulation output. 
# Specifically, checking to see that the total number of unique species in all.populations
# is the same as the number of species in phylo.out.
# The output contains a row for each sim with the two measures of richness and their difference.

file_dir = 'C:/SENCoutput'


phylo.allpops.check = function(file_dir) {
  require(ape)
  
  files = list.files(file_dir)
  allpops.files = files[grep('all.pops',files)]
  phylo.files = files[grep('phylo',files)]
  
  sims = as.vector( sapply(allpops.files, function(x)
                      as.numeric(strsplit(strsplit(x, "sim")[[1]][2] , "\\.csv")[[1]][1])) )
  
  check = data.frame(sim = NA, allpops.S = NA, phylo.S = NA)
  for (i in sims) {
    pops = read.csv(paste(file_dir,'/SENC_all.pops_sim',i,'.csv',sep=''), header=T)
    phy = read.tree(paste(file_dir,'/SENC_phylo_sim',i,'.tre',sep=''))
    if (nrow(pops) > 0) {
      check = rbind(check, c(i, length(unique(pops$spp.name)), length(phy$tip.label)))
    } else {
      check = rbind(check, c(i, NA, NA))
    }
  }
  check = check[-1,]
  check$diff = check$phylo.S - check$allpops.S
  check = check[order(check$diff,decreasing = T),]
  return(check)
}
