
#Function for unzipping output files, and storing them in a list
# --output_dir is the directory where the zipped output files reside

output.unzip = function(output_dir, sim_ID) {
  require(ape)
  
  #check that zip file exists (returns value of 0 when mode=0 if file exists)
  if( file.access(paste(output_dir,'/senc.out.',sim_ID,'.zip',sep=''), mode=0) == 0 ) {
    unzipped.filenames = unzip(paste(output_dir,'/senc.out.',sim_ID,'.zip',sep=''), exdir = output_dir)
    all.populations = read.csv(unzipped.filenames[1], header=T)
    params.out = read.csv(unzipped.filenames[2],header=T)
    phylo.out = read.tree(unzipped.filenames[3])
    time.richness = read.csv(unzipped.filenames[4], header=T)
    sim.results = list(all.populations=all.populations, 
                       params.out=params.out, 
                       phylo.out=phylo.out, 
                       time.richness=time.richness)
  } else if (file.access(paste(output_dir, '/SENC_all.pops_sim', sim_ID, '.csv', sep = ''), mode = 0) == 0) { 
    #if sim output data are not zipped, then just read in raw files
    all.populations = read.csv(paste(output_dir, '/SENC_all.pops_sim', sim_ID, '.csv', sep = ''), header=T)
    params.out = read.csv(paste(output_dir, '/SENC_params.out_sim', sim_ID, '.csv', sep = ''), header=T)
    phylo.out = read.tree(paste(output_dir, '/SENC_phylo_sim', sim_ID, '.tre', sep = ''))
    time.richness = read.csv(paste(output_dir, '/SENC_time.rich_sim', sim_ID, '.csv', sep = ''), header=T)
    sim.results = list(all.populations=all.populations, 
                       params.out=params.out, 
                       phylo.out=phylo.out, 
                       time.richness=time.richness)
  } else {
    print("No simulation output exists for that sim ID")
    sim.results = NULL
  }
    return(sim.results)
}

