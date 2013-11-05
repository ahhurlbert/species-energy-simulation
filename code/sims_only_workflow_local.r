# Simulation workflow
user="Allen"
partial.analysis=1


# (1) set wd for location to write results
if (user=="Allen") {
  setwd("C:/SENCoutput")
}
if (user=="James") {
  setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204")
}


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
library(parallel);
library(doParallel);


source('code/senc_sim_fun.r')
source('code/make.phylo.jimmy.fun.r')

#(3) read in master simulation matrix with chosen parameter combinations
# then add fields for storing output summary
sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv", header=T)
sim.matrix$n.regions = NA
sim.matrix$extant.S = NA
sim.matrix$extinct.S = NA
sim.matrix$skipped.clades = NA
sim.matrix$skipped.times = NA
sim.matrix$BK.reg = NA
sim.matrix$BK.env = NA

package.vector = c('ape','permute','nlme','vegan','picante','mvtnorm','caper','paleotree','plyr','phytools');

cl = makeCluster(2);
registerDoParallel(cl);


#(4) start analyses based on value of 'sim' which draws parameter values from sim.matrix
if (partial.analysis == 0) {which.sims = 1:max(sim.matrix$sim.id)};
if (partial.analysis == 1) {which.sims = 3125:3324}; # which.sims = c(read.csv(paste(analysis_dir,"/sims.to.analyze.csv",sep=""))$x)

foo = foreach(sim=which.sims,.packages = package.vector,.combine='rbind') %dopar% {
  
  #(4) call simulation
  sim.start = date(); print(sim.start);
  sim.results = senc_sim_fun(sim.matrix=sim.matrix,sim=sim)
  sim.end = date(); print(sim.end);
  write.csv(sim.results$all.populations,paste('SENC_all.pops_sim',sim,'.csv',sep=''), quote=F,row.names=F)
  write.csv(sim.results$time.richness,paste("SENC_time.rich_sim",sim,".csv",sep=""),quote=F,row.names=F);
  write.tree(sim.results$phylo.out,paste("SENC_phylo_sim",sim,".tre",sep=""));
  write.csv(sim.results$params.out,paste("SENC_params.out_sim",sim,".csv",sep=""),quote=F,row.names=F);
  
}
