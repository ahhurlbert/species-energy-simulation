# R script to be called from the shell that will run the specified simulations
# Typical call on the UNC KillDevil cluster (which uses LSF) would be
#
#   $bsub -o out.%J -n 100 -a openmpi mpirun Rscript run_sim_on_cluster.r 3765 3864 30
#
# which would use 100 nodes to run sims 3765:3864 at 30 evenly spaced time points

source('code/run_sim.r')
source('code/analyze_sim.r')
source('code/senc_analysis_fun.r')
source('code/supplemental_analysis_functions.r')

args = commandArgs(trailingOnly = TRUE)

# Specify a vector of sim.id's to run in the command line call
# --if a simple sequence, then arguments are simply the first and second sims
# --if more complex set, then arguments are paired to represent subsequences, and then catenated
# --e.g. R script run_sim_on_cluster.r 1 18 26 29 5 would run c(1:18, 26:29) at 5 points in time
which_sims = c()
for (i in 1:((length(args) - 1)/2)) {
  which_sims = c(which_sims, as.numeric( args[2*i - 1] ) : as.numeric( args[2*i]) )
}

times = as.numeric(args[length(args)])

run_sim(which_sims, "SENC_Master_Simulation_Matrix.csv", local = F)

analyze_sim(which_sims, "SENC_Master_Simulation_Matrix.csv", local = F, sim_dir = "new", num.of.time.slices = times)
