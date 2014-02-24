# Use this function to analyze simulations by specifying:
#   which_sims:  the simulation IDs to run
#   sim_matrix_filename:  a .csv file specifying parameter combinations for each sim ID
#   local:       TRUE if running on a local machine, FALSE if running on a remote cluster
#   num_cores:   number of processors available for parallel processing
#   root.only:   analyze root clade only (1) or all subclades (0)
#   num.of.time.slices:  number of points in time to analyze
#   which.time.slices:   which points in time to analyze if irregularly spaced
#   time.sequence:       which points in time to analyze if regularly spaced
#   min.num.spp:         min num of spp in a clade needed to proceed with analysis
#   sim_dir:             directory where raw simulation output reside ('archived' for example
#                          output from the published paper, 'new' for simulation output created
#                          by the user)

analyze_sim = function(which_sims, 
                       sim_matrix_filename = "SENC_Master_Simulation_Matrix.csv", 
                       local = T, 
                       num_cores = 1, 
                       root.only = 1, 
                       num.of.time.slices = 1, 
                       which.time.slices = NA, 
                       time.sequence = NA,  
                       min.num.spp = 8,
                       sim_dir = 'archived') 
{

  require(ape)
  require(caper)
  require(paleotree)
  require(plyr)
  require(foreach)
  
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
  writeLines(c(""), "analysis_output/analysis_log.txt")
  sink("analysis_output/analysis_log.txt", append = T)
  
  package.vector = c('ape','caper','paleotree','plyr')
  export.vector = c("analysis", "output.unzip", "regional.calc",
                    "xregion.analysis", "extinct.calc", "maxlik.betasplit.AH")
  
  foo = foreach(sim = which_sims, .packages = package.vector, .combine = "rbind",
                .export = export.vector) %dopar% {
                                
    sink("analysis_output/analysis_log.txt", append = T)
    cat(paste("Starting sim", sim, ",", date(), "\n"))
                                
    analysis(sim, sim_matrix, root.only = root.only, num.of.time.slices = num.of.time.slices, 
             which.time.slices = which.time.slices, time.sequence = time.sequence, 
             min.num.spp = min.num.spp, sim_dir = sim_dir)
                                
  } #end foreach

  #stop or close clusters
  if(local) {
    stopCluster(cl)
  } else {
    closeCluster(cl)
    mpi.finalize()
  }
  
} #end function
