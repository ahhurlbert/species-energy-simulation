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

source('reg_calc_and_analysis.r')
source('senc_sim_fun.r');
source('make.phylo.jimmy.fun.r');

#(3) read in master simulation matrix with chosen parameter combinations
sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));

#(4) start simulation based on value of 'sim' which draws parameter values from sim.matrix
output = numeric();
num.of.time.slices = 5;
      
#call simulation, with output being the equivalent of the all.populations dataframe

sim.start = date(); print(sim.start);
sim.results = senc_sim_fun(sim.matrix=sim.matrix,sim=sim)
sim.end = date(); print(sim.end);
    
timeslices = as.integer(round(seq(max(sim.results$time.richness$time)/num.of.time.slices,max(sim.results$time.richness$time),length=num.of.time.slices),digits=0));
max.time.actual = max(sim.results$time.richness$time);

for (t in timeslices) {
	sub.species = subset(sim.results$all.populations,time.of.sp.origin < t & time.of.sp.extinction > t);
	tips.to.drop = as.character(sim.results$phylo.out$tip.label[which(is.element(sim.results$phylo.out$tip.label,as.character(sub.species$spp.name))==F)]);
  	sub.phylo = drop.tip(sim.results$phylo.out,tips.to.drop);
	temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1]); temp.root.time;
	most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1])))]; most.recent.spp;
	extinct.time.most.recent = unique(sim.results$all.populations$time.of.sp.extinction[sim.results$all.populations$spp.name==most.recent.spp]); extinct.time.most.recent;
	sub.phylo$root.time = temp.root.time + max(c(0,max.time.actual - extinct.time.most.recent)); sub.phylo$root.time;
	sub.phylo = collapse.singles(timeSliceTree(sub.phylo,sliceTime=(max.time.actual - t),plot=F,drop.extinct = T));
	num.of.spp = length(sub.phylo$tip.label);
	rm('tips.to.drop','temp.root.time');

	for (c in (num.of.spp+1):max(sub.phylo$edge)) {

		#pull out list of species names belonging to each subclade
 		sub.clade = clade.members(c, sub.phylo, tip.labels=T)
 		subset.populations = subset(sim.results$all.populations, spp.name %in% as.numeric(sub.clade));
			
		#sub.populations is the subset of populations specific to a particular clade and timeslice
		sub.populations = subset(subset.populations, time.of.origin < t & time.of.extinction > t)

		if (length(unique(sub.populations$spp.name))>1) {

			#sub.clade.phylo is a specific simulation clade pulled from the phylogeny that was sliced at timeslice t
			tips.to.drop = as.character(sub.phylo$tip.label[which(is.element(sub.phylo$tip.label,as.character(sub.populations$spp.name))==F)]);
  			sub.clade.phylo = drop.tip(sub.phylo,tips.to.drop);
			sub.clade.phylo$root.time = max(dist.nodes(sub.clade.phylo)[1:Ntip(sub.clade.phylo),Ntip(sub.clade.phylo) + 1]); sub.clade.phylo$root.time;
			sub.clade.phylo$origin.time = t - sub.clade.phylo$root.time; sub.clade.phylo$origin.time;

			if (identical(sort(as.integer(unique(sub.populations$spp.name))) , sort(as.integer(sub.clade.phylo$tip.label)))==F ) {print(c(c,t,'Error: trimmed phylogeny does not contain the correct species')); break} else{}; 

    			reg.summary = regional.calc(sub.populations[,c('region','spp.name','time.of.origin','reg.env')], sub.clade.phylo, t);

    			corr.results = xregion.analysis(reg.summary)
          
  			output = rbind(output, cbind(sim=sim,clade.id = c, time = t, corr.results))
			print(c(c,t,date(),length(sub.clade.phylo$tip.label),extinct.time.most.recent));

 		} else{};
	} # end sub clade for loop
        
}; # end timeslice loop
      
#write all of this output to files
write.csv(output,paste("SENC_Stats_sim",sim,".csv",sep=""),quote=F,row.names=F);
analysis.end = date();
print(c(warnings(),sim.start,sim.end,analysis.end));


####################################################
# Simulation summary plots
####################################################
lat.grad.time.plot(sim.results, numslices = 10)
clade.origin.corr.plot(output, sim.results$params.out)
clade.exmpl.figs(sim.results, clade.slices=6, seed=0)
