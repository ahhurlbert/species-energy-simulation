# Use this function to run simulations by specifying:
#   which_sims:  the simulation IDs to run
#   sim_matrix:  a .csv file specifying parameter combinations for each sim ID
#   local:       TRUE if running on a local machine, FALSE if running on a remote cluster
#   num_cores:   number of processors available for parallel processing

source("code/senc_sim_fun.r")

run_sim = function(which_sims, 
                   sim_matrix_filename = "SENC_Master_Simulation_Matrix.csv", 
                   local = T, 
                   num_cores = 1) 
{

  require(ape)
  require(parallel)
  require(doParallel)
  
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  
  sim_matrix = read.csv(sim_matrix_filename, header = T)

  if(local) {
    foo = foreach(sim = which_sims, .packages = 'ape', .combine = "rbind", .export = "senc_sim_fun") %dopar% {
      
      sim_start = date()
      print(sim_start)
      sim.results = senc_sim_fun(sim_matrix, sim = sim)
      sim_end = date()
      print(sim_end)
      
    } #end foreach
    
  } #end if(local)

  if(!local) {
    #Pass in arguments from shell command line? James, is something like this necessary?
  } #end if(!local)
  
} #end function

