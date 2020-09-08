test_fun_mpi = function(which_sims) 
{
    require(ape)
    require(parallel)
    require(doParallel)
    require(doMPI)
        
        # create and register a doMPI cluster if necessary
        if (!identical(getDoParName(), 'doMPI')) {
          cl <- startMPIcluster()
          registerDoMPI(cl)
    }
  
  #Create a log file for checking simulation progress
  writeLines(c(""), "raw_sim_output/test_log.txt")
    
  foo = foreach(sim = which_sims, .packages = 'ape', .combine = "rbind")  %dopar% {
      
    sink("raw_sim_output/log.txt", append = T)
    cat(paste("Starting sim", sim, ",", date(), "\n"))
    sim.results = sim
    print(sim)
    cat(paste("Finished sim", sim, ",", date(), "\n"))
    sink()
  } #end foreach
    
    closeCluster(cl)
    mpi.finalize()
  
} #end function
