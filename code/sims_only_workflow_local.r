# Simulation workflow

#(1) load libraries and simulation and phylogeny-making functions
library(ape)

# Info for parallelizing processing specify number of processors in makeCluster
package.vector = c('ape')
cl = makeCluster(2)
registerDoParallel(cl)

source("code/senc_sim_fun.r")
source("code/make.phylo.jimmy.fun.r")


#(2) read in master simulation matrix with chosen parameter combinations
sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv", header=T)


#(3) Specify which simulation IDs to run, as given in the sim.matrix
which.sims = 3125:3125

foo = foreach(sim = which.sims, .packages = package.vector, .combine = "rbind") %dopar% {
  
  #(4) call simulation, output is written to the raw_sim_output folder
  sim.start = date()
  print(sim.start)
  sim.results = senc_sim_fun(sim.matrix = sim.matrix, sim = sim)
  sim.end = date()
  print(sim.end)
  write.csv(sim.results$all.populations,paste("raw_sim_output/SENC_all.pops_sim", sim, ".csv", sep = ""), quote = F, row.names = F)
  write.csv(sim.results$time.richness,paste("raw_sim_output/SENC_time.rich_sim", sim, ".csv", sep = ""), quote = F, row.names = F)
  write.tree(sim.results$phylo.out,paste("raw_sim_output/SENC_phylo_sim", sim, ".tre", sep = ""))
  write.csv(sim.results$params.out,paste("raw_sim_output/SENC_params.out_sim", sim, ".csv", sep = ""), quote = F, row.names = F)
  
}
