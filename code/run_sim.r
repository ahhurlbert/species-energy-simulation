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

    #Create a log file for checking simulation progress
    writeLines(c(""), "log.txt")
    counter = 0
    
    foo = foreach(sim = which_sims, .packages = 'ape', .combine = "rbind", .export = "senc_sim_fun") %dopar% {
      
      counter = counter + 1
      sink("raw_sim_output/log.txt", append = T)
      cat(paste("Starting sim", sim, "(", counter, "out of", length(which_sims), "),", date(), "\n"))
      sim.results = senc_sim_fun(sim_matrix, sim = sim)
      
    } #end foreach
    
  } #end if(local)

  if(!local) {
    #Pass in arguments from shell command line? James, is something like this necessary?
  } #end if(!local)
  
} #end function

