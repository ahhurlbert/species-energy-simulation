# R script to be called from the shell that will run the specified simulations

# Test running on cluster, but not in parallel

source('code/run_sim.r')

args = commandArgs(trailingOnly = TRUE)

# Specify a vector of sim.id's to run in the command line call
# --if a simple sequence, then arguments are simply the first and last sims
# --if more complex set, then arguments are paired to represent subsequences, and then catenated
# --e.g. R script run_sim_on_cluster.r 1 18 26 29  would run c(1:18, 26:29)
which_sims = c()
for (i in 1:(length(args)/2)) {
  which_sims = c(which_sims, as.numeric( args[2*i - 1] ) : as.numeric( args[2*i]) )
}

run_sim(which_sims, "test_matrix.csv", local = T)