# Simulation workflow

# (1) set wd for location to write results

setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130111")

#(2) load libraries simulation and analysis functions
library(ape);
library(permute);
library(nlme);
library(vegan);
library(picante);
library(mvtnorm);
library(caper);
library(paleotree);
library(plyr);

source('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation/senc_sim_fun.r');
source('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation/make.phylo.jimmy.fun.r');

#(3) read in master simulation matrix with chosen parameter combinations
sim.matrix = as.data.frame(read.csv("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation/SENC_Master_Simulation_Matrix.csv",header=T));

#(4) call simulation
sim = 2604
sim.start = date(); print(sim.start);
sim.results = senc_sim_fun(sim.matrix=sim.matrix,sim=sim)
sim.end = date(); print(sim.end);