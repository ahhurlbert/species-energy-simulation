# Code for analyzing the simulated phylogenies from 
# Hurlbert & Stegen 2014a,b that were generated under
# different eco-evolutionary dynamics using Morlon et al.'s
# (2010, PLoS Biology) RPANDA package

# Set working directory to species-energy-simulation repo

library(RPANDA)

# Read in trees
sims = 4665

output = c()
for (i in sims) {
  treename = paste('raw_sim_output/sim', i, '_out/SENC_phylo_sim', i, '.tre', sep = '')
  tree = read.tree(treename)
  out = fit_coal_cst(tree, tau0=1.e-3, gamma=-1, cst.rate=FALSE)
}