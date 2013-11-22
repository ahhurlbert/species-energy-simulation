# load simulation and analysis functions
library(mnormt)
library(rgl)
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
library(apTreeshape)
library(foreach)
library(doParallel)

package.vector = c('ape','permute','nlme','vegan','picante','mvtnorm','caper','paleotree','plyr','phytools','apTreeshape')

source('code/reg_calc_and_analysis.r')
source('code/make.phylo.jimmy.fun.r')
source('code/lat.grad.time.plot.r')
source('code/clade.origin.corr.plot.r')
source('code/clade.exmpl.figs.r')
source('code/extinct.calc.r')
source('code/unzipping_files.r')
source('code/maxlik.betasplit.AH.r')
source('code/analysis_workflow_Olympus.r')

# read in master simulation matrix with chosen parameter combinations;
# then add fields for storing output summary
sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv",header=T)

#####################################################################
# Specify sims to run

#sim,  simulation ID to analyze
#sim.matrix,  matrix with simulation parameters for all simulations
#local,   analysis on a local machine (1) or cluster (0)
#root.only, analyze root clade only (1) or all subclades (1)
#num.of.time.slices,  number of points in time to analyze (if 1, then the end of the simulation)
#which.time.slices,  which points in time to analyze if irregularly spaced (a vector)
#time.sequence,   which points in time to analyze if regularly spaced (a vector)
#min.num.spp,   minimum number of species in a clade needed to proceed with analysis
#num.processors,  number of processors for parallel processing on local machine
  
analysis(sim = ,
         sim.matrix,
         local = 1,
         root.only = 1,
         num.of.time.slices = 1,
         which.time.slices = NA,
         time.sequence = NA,
         min.num.spp = 8,
         num.processors = 2)
