# R script to be called from the shell that will run the specified simulations
# Typical call on the UNC KillDevil cluster (which uses LSF) would be
#
#   $bsub -o out.%J -n 100 -a openmpi mpirun Rscript run_sim_on_cluster.r 3765 3864
#
# which would use 100 nodes to run sims 3765:3864 

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

run_sim(which_sims, "SENC_Master_Simulation_Matrix.csv", local = FALSE)
