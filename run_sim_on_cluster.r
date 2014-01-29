library(ape)
library(permute)
library(nlme)
library(vegan)
library(picante)
library(mvtnorm)
library(caper)
library(paleotree)
library(plyr)
library(phytools)
library(foreach)
library(doParallel)
library(apTreeshape)
library(mnormt)
library(rgl)
library(apTreeshape)
library(R.utils)
library(ncdf)


source('code/run_sim.r')
source('code/analyze_sim.r')
source('code/senc_sim_fun.r')
source('code/senc_analysis_fun.r')
source('code/supplemental_analysis_functions.r')

args = commandArgs(trailingOnly = TRUE)

which_sims = c()
for (i in 1:(length(args)/2)) {
  which_sims = c(which_sims, as.numeric( args[2*i - 1] ) : as.numeric( args[2*i]) )
}

run_sim(which_sims, "SENC_Master_Simulation_Matrix.csv", local = F)