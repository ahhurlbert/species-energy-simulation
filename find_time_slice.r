#!/usr/bin/env Rscript

#sim = commandArgs();
#sim = as.numeric(sim[length(sim)]);

# Choose number of time slices per simulation to analyze
num.of.time.slices = 1;
# Set minimum number of species in a clade needed to proceed with analysis
min.num.spp = 8;

Allen = 0;
partial.analysis = 1; # toggle to determine whether we're looking at all sims or just some

#New parameter for taking into account which of us is running this code
if(Allen==1) {
  setwd('c:/documents and settings/hurlbert/species-energy-simulation')
  Rlib.location = "C:/program files/R/R-2.15.2/library"
  sim_dir = "C:/SENCoutput"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses"
} else {
  setwd('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation')
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever all of your zipped output files are
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever you want to store the results of these analyses
}

# Simulation workflow

#(2) load simulation and analysis functions
if (Allen==1){
  library(ape,lib.loc=Rlib.location);
  library(permute,lib.loc=Rlib.location);
  library(nlme,lib.loc=Rlib.location);
  library(vegan,lib.loc=Rlib.location);
  library(picante,lib.loc=Rlib.location);
  library(mvtnorm,lib.loc=Rlib.location);
  library(caper,lib.loc=Rlib.location);
  library(paleotree,lib.loc=Rlib.location);
  library(plyr,lib.loc=Rlib.location);
  library(phytools, lib.loc=Rlib.location);
  library(foreach,lib.loc=Rlib.location);
  library(doParallel,lib.loc=Rlib.location);
} else {
  library(ape);
  library(permute);
  library(nlme);
  library(vegan);
  library(picante);
  library(mvtnorm);
  library(caper);
  library(paleotree);
  library(plyr);
  library(phytools);
  library(foreach);
  library(doParallel);
}

source('reg_calc_and_analysis.r');
source('make.phylo.jimmy.fun.r');
source('lat.grad.time.plot.r');
source('clade.origin.corr.plot.r');
source('clade.exmpl.figs.r');
source('extinct.calc.r');
source('unzipping_files.r');

#(3) read in master simulation matrix with chosen parameter combinations;
# then add fields for storing output summary
sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));
head(sim.matrix);

tropical.sims = sim.matrix$sim.id[sim.matrix$reg.of.origin == 'tropical' & sim.matrix$carry.cap == 'on' & sim.matrix$energy.gradient == 'on' & sim.matrix$sim.id > 3464]; length(tropical.sims);
temperate.sims = sim.matrix$sim.id[sim.matrix$reg.of.origin == 'temperate' & sim.matrix$carry.cap == 'on' & sim.matrix$energy.gradient == 'on' & sim.matrix$sim.id > 3464]; length(temperate.sims);

trop.orig.extreme = 3;
temp.orig.extreme = 8;

trop.times = numeric();

for (i in 1:length(tropical.sims)) {
  
  sim = tropical.sims[i];
  sim.results = output.unzip(sim_dir,sim)
  time.richness = sim.results$time.richness
  trop.times = c(trop.times,min(time.richness$time[time.richness$region == trop.orig.extreme & time.richness$spp.rich >= 5]));
  print(i)
  
}

hist(trop.times); quantile(trop.times);

temp.times = numeric();

for (i in 1:length(temperate.sims)) {
  
  sim = temperate.sims[i];
  sim.results = output.unzip(sim_dir,sim)
  time.richness = sim.results$time.richness
  temp.times = c(temp.times,min(time.richness$time[time.richness$region == temp.orig.extreme & time.richness$spp.rich >= 5]));
  print(i)
  
}

hist(temp.times); quantile(temp.times);
#0%    25%    50%    75%   100% 
#2508.0 4442.5 5459.0 6359.5 8674.0 



