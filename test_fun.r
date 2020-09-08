test_fun = function(which_sims, local = T, num_cores = 1) 
{
    require(ape)
    require(parallel)
    require(doParallel)
    cl = makeCluster(num_cores)
    registerDoParallel(cl)
  
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
    
    stopCluster(cl)
  
} #end function
