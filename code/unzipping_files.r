
#Function for unzipping output files, and storing them in a list
# --output_dir is the directory where the zipped output files reside

output.unzip = function(raw_sim_output, sim_ID) {
  require(ape)
  
  
  #check that zip file exists (returns value of 0 when mode=0 if file exists)
  if( file.access(paste(raw_sim_output,'/senc.out.',sim_ID,'.zip',sep=''), mode=0) == 0 ) {
    unzipped.filenames = unzip(paste(raw_sim_output,'/senc.out.',sim_ID,'.zip',sep=''), exdir = raw_sim_output)
    all.populations = read.csv(unzipped.filenames[1], header=T)
    params.out = read.csv(unzipped.filenames[2],header=T)
    phylo.out = read.tree(unzipped.filenames[3])
    time.richness = read.csv(unzipped.filenames[4], header=T)
  } else { #if sim output data are not zipped, then just read in raw files
    all.populations = read.csv(paste(raw_sim_output, '/SENC_all.pops_sim', sim_ID, '.csv', sep = ''), header=T)
    params.out = read.csv(paste(raw_sim_output, '/SENC_params.out_sim', sim_ID, '.csv', sep = ''), header=T)
    phylo.out = read.tree(paste(raw_sim_output, '/SENC_phylo_sim', sim_ID, '.tre', sep = ''))
    time.richness = read.csv(paste(raw_sim_output, '/SENC_time.rich_sim', sim_ID, '.csv', sep = ''), header=T)
  }
  sim.results = list(all.populations=all.populations, 
                     params.out=params.out, 
                     phylo.out=phylo.out, 
                     time.richness=time.richness)
  return(sim.results)
}

