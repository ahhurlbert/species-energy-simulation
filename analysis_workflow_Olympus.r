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

source('reg_calc_and_analysis.r');
source('make.phylo.jimmy.fun.r');
source('lat.grad.time.plot.r');
source('clade.origin.corr.plot.r');
source('clade.exmpl.figs.r');
source('extinct.calc.r');

#(3) read in master simulation matrix with chosen parameter combinations
sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));

#(4) start analyses based on value of 'sim' which draws parameter values from sim.matrix
output = numeric();
num.of.time.slices = 5;
      
# (5) read in simulation results for specified simulation

all.populations = read.csv(paste("SENC_all.pops_sim",sim,".csv",sep=""),header=T);
time.richness = read.csv(paste("SENC_time.rich_sim",sim,".csv",sep=""),header=T);
phylo.out = read.tree(paste("SENC_phylo_sim",sim,".tre",sep=""));
params.out = read.csv(paste("SENC_params.out_sim",sim,".csv",sep=""),header=T);

sim.results = list(all.populations=all.populations,time.richness=time.richness,phylo.out=phylo.out,params.out=params.out);

timeslices = as.integer(round(seq(max(sim.results$time.richness$time)/num.of.time.slices,max(sim.results$time.richness$time),length=num.of.time.slices),digits=0));
max.time.actual = max(sim.results$time.richness$time);

for (t in timeslices) {
	sub.species = as.character(unique(subset(all.populations,time.of.sp.origin < t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);

  #Some species may be extant globally but in our boundary regions (0,11) only; 
  #we need to eliminate species that are not extant within regions 1-10
  extant.ornot = aggregate(all.populations$extant,by=list(all.populations$spp.name),sum)
  extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
  sub.species2 = sub.species[!sub.species %in% extinct.species]
	tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
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

		if (length(unique(sub.populations$spp.name))>=10) {

			#sub.clade.phylo is a specific simulation clade pulled from the phylogeny that was sliced at timeslice t
			tips.to.drop = as.character(sub.phylo$tip.label[which(is.element(sub.phylo$tip.label,as.character(sub.populations$spp.name))==F)]);
  		sub.clade.phylo = drop.tip(sub.phylo,tips.to.drop);
			sub.clade.phylo$root.time = max(dist.nodes(sub.clade.phylo)[1:Ntip(sub.clade.phylo),Ntip(sub.clade.phylo) + 1]); sub.clade.phylo$root.time;
			sub.clade.phylo$origin.time = t - sub.clade.phylo$root.time; sub.clade.phylo$origin.time;

			if (identical(sort(as.integer(unique(sub.populations$spp.name))) , sort(as.integer(sub.clade.phylo$tip.label)))==F ) {print(c(c,t,'Error: trimmed phylogeny does not contain the correct species')); break} else{}; 

    	reg.summary = regional.calc(sub.populations[,c('region','spp.name','time.of.origin','reg.env','extant')], sub.clade.phylo, t);
	
      #Note that extinction calculation must be done on subset.populations, not sub.populations
      extinction = extinct.calc(subset.populations, timeslice=t)
      reg.summary2 = merge(reg.summary,extinction[,c('region','extinction.rate')],by='region')
      
    	corr.results = xregion.analysis(reg.summary2)
          
  		#Pybus & Harvey (2000)'s gamma statistic
      Gamma.stat = gammaStat(sub.clade.phylo)
      
      output = rbind(output, cbind(sim=sim,clade.id = c, time = t, corr.results, gamma.stat = Gamma.stat))
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
clade.exmpl.figs(sim.results, output, clade.slices=6, seed=0)
