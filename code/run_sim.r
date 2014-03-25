# Use this function to run simulations by specifying:
#   which_sims:  the simulation IDs to run
#   sim_matrix_filename:  a .csv file specifying parameter combinations for each sim ID
#   local:       TRUE if running on a local machine, FALSE if running on a remote cluster
#   num_cores:   number of processors available for parallel processing

run_sim = function(which_sims, 
                   sim_matrix_filename = "SENC_Master_Simulation_Matrix.csv", 
                   local = T, 
                   num_cores = 1) 
{

  require(ape)
  source('code/senc_sim_fun.r')
  source('code/supplemental_analysis_functions.r')
  
  if(local) {
    require(parallel)
    require(doParallel)
    cl = makeCluster(num_cores)
    registerDoParallel(cl)
  } else {
    require(doMPI)
    
    # create and register a doMPI cluster if necessary
    if (!identical(getDoParName(), 'doMPI')) {
      cl <- startMPIcluster()
      registerDoMPI(cl)
    }
  }
    
  sim_matrix = read.csv(sim_matrix_filename, header = T)
  
  #Create a log file for checking simulation progress
  writeLines(c(""), "raw_sim_output/log.txt")
    
  foo = foreach(sim = which_sims, .packages = 'ape', .combine = "rbind", .export =
                  "senc_sim_fun")  %dopar% {
      
    sink("raw_sim_output/log.txt", append = T)
    cat(paste("Starting sim", sim, ",", date(), "\n"))
    sim.results = senc_sim_fun(sim_matrix, sim = sim)
    cat(paste("Finished sim", sim, ",", date(), "\n"))  
  } #end foreach
    
  if(local) {
    stopCluster(cl)
  } else {
    closeCluster(cl)
    mpi.finalize()
  }

} #end function

