#!/usr/bin/env Rscript

sim = commandArgs();
sim = as.numeric(sim[length(sim)]);

# Simulation workflow

#(2) load simulation and analysis functions
library(ape,lib.loc="/pic/people/steg815/Rlibs");
library(permute,lib.loc="/pic/people/steg815/Rlibs");
library(nlme,lib.loc="/pic/people/steg815/Rlibs");
library(vegan,lib.loc="/pic/people/steg815/Rlibs");
library(picante,lib.loc="/pic/people/steg815/Rlibs");
library(mvtnorm,lib.loc="/pic/people/steg815/Rlibs");
library(caper,lib.loc="/pic/people/steg815/Rlibs");
library(paleotree,lib.loc="/pic/people/steg815/Rlibs");
library(plyr,lib.loc="/pic/people/steg815/Rlibs");

source('senc_sim_fun.r');
source('make.phylo.jimmy.fun.r');

#(3) read in master simulation matrix with chosen parameter combinations
sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));

#call simulation, with output being the equivalent of the all.populations dataframe

sim.start = date(); print(sim.start);
sim.results = senc_sim_fun(sim.matrix=sim.matrix,sim=sim)
sim.end = date(); print(sim.end);