
#Function for unzipping output files, and storing them in a list
# --output_dir is the directory where the zipped output files reside

output.unzip = function(sim_dir, sim_ID) {
  require(ape)
  
  
  #check that file exists (returns value of 0 when mode=0 if file exists)
  if( file.access(paste(sim_dir,'/senc.out.',sim_ID,'.zip',sep=''), mode=0) == 0 ) {
    unzipped.filenames = unzip(paste(sim_dir,'/senc.out.',sim_ID,'.zip',sep=''), exdir = sim_dir)
    all.populations = read.csv(unzipped.filenames[1], header=T)
    params.out = read.csv(unzipped.filenames[2],header=T)
    phylo.out = read.tree(unzipped.filenames[3])
    time.richness = read.csv(unzipped.filenames[4], header=T)
    sim.results = list(all.populations=all.populations, 
                       params.out=params.out, 
                       phylo.out=phylo.out, 
                       time.richness=time.richness)
    return(sim.results)
  } else { 
    return(NULL)}
}

