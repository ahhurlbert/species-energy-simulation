
#Function for unzipping output files
# --output_dir is the directory where the zipped output files reside

output.unzip = function(output_dir, sim_ID) {
  unzipped.filenames = unzip(paste(output_dir,'/senc.out.',sim_ID,'.zip',sep=''))
  all.populations = read.csv(unzipped.filenames[1], header=T)
  params.out = read.csv(unzipped.filenames[2],header=T)
  phylo.out = read.tree(unzipped.filenames[3])
  time.richness = read.csv(unzipped.filenames[4], header=T)
  sim.results = list(all.populations=all.populations, 
                     params.out=params.out, 
                     phylo.out=phylo.out, 
                     time.richness=time.richness)
  return(sim.results)
}

