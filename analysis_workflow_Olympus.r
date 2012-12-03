#!/usr/bin/env Rscript

#sim = commandArgs();
#sim = as.numeric(sim[length(sim)]);

#New parameter for taking into account which of us is running this code
if(Allen==1) {
  setwd('c:/documents and settings/hurlbert/species-energy-simulation')
  Rlib.location = "C:/program files/R/R-2.15.1/library"  
  sim_dir = "C:/SENCoutput"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses"
}
else {
  Rlib.location = "/pic/people/steg815/Rlibs"            
  sim_dir = ............ #wherever all of your zipped output files are
  analysis_dir = .......... #wherever you want to store the results of these analyses
}

# Simulation workflow

#(2) load simulation and analysis functions
library(ape,lib.loc=Rlib.location);
library(permute,lib.loc=Rlib.location);
library(nlme,lib.loc=Rlib.location);
library(vegan,lib.loc=Rlib.location);
library(picante,lib.loc=Rlib.location);
library(mvtnorm,lib.loc=Rlib.location);
library(caper,lib.loc=Rlib.location);
library(paleotree,lib.loc=Rlib.location);
library(plyr,lib.loc=Rlib.location);

source('reg_calc_and_analysis.r');
source('make.phylo.jimmy.fun.r');
source('lat.grad.time.plot.r');
source('clade.origin.corr.plot.r');
source('clade.exmpl.figs.r');
source('extinct.calc.r');
source('unzipping_files.r');

#(3) read in master simulation matrix with chosen parameter combinations, and define number of 
#    time slices to analyze
sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));
num.of.time.slices = 5;


#(4) start analyses based on value of 'sim' which draws parameter values from sim.matrix
which.sims = 1:max(sim.matrix$sim.id)

for (sim in which.sims) {
  rm(list=c('all.populations', 'time.richness', 'phylo.out', 'params.out', 'output', 'sim.results'))
  output = numeric();

  # (5) read in simulation results for specified simulation from the output zip file
  sim.results = output.unzip(sim_dir,sim)
  all.populations = sim.results$all.populations
  time.richness = sim.results$time.richness
  phylo.out = sim.results$phylo.out
  params.out = sim.results$params.out

  
  timeslices = as.integer(round(seq(max(time.richness$time)/num.of.time.slices,max(time.richness$time),length=num.of.time.slices),digits=0));
  max.time.actual = max(time.richness$time);
  
  for (t in timeslices) {
    # vector of species in existence at time t
    sub.species = as.character(unique(subset(all.populations,time.of.sp.origin < t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);
    
    # Some species may be extant globally (extant==1) but in our boundary regions (0,11) only; 
    # we need to eliminate species that are not extant within regions 1-10 (which is all that is 
    # reflected in the all.populations dataframe)
    # FIXME: 
    # Add more explanatory comments justifying why we don't need to consider species that existed
    # at time t but went extinct before the present.
    # In some cases (e.g. sim 2, t=6000), tips.to.drop includes all tips and so sub.phylo is empty.
    # Does it make sense for this to ever happen? If not, fix it. 
    # If so, need to provide an if-else error catch both in the creation of sub.phylo, 
    # and of sub.clade.phylo inside the clade loop. (Sim 3, t = 156 bonks at that point)
    extant.ornot = aggregate(all.populations$extant,by=list(all.populations$spp.name),sum)
    extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
    sub.species2 = sub.species[!sub.species %in% extinct.species]
    tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
    sub.phylo = drop.tip(phylo.out,tips.to.drop);
    
    temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1]); temp.root.time;
    most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1])))]; most.recent.spp;
    extinct.time.most.recent = unique(all.populations$time.of.sp.extinction[all.populations$spp.name==most.recent.spp]); extinct.time.most.recent;
    sub.phylo$root.time = temp.root.time + max(c(0,max.time.actual - extinct.time.most.recent)); sub.phylo$root.time;
    sub.phylo = collapse.singles(timeSliceTree(sub.phylo,sliceTime=(max.time.actual - t),plot=F,drop.extinct = T));
    num.of.spp = length(sub.phylo$tip.label);
    rm('tips.to.drop','temp.root.time');
    
    for (c in (num.of.spp+1):max(sub.phylo$edge)) {
      
      #pull out list of species names belonging to each subclade
      sub.clade = clade.members(c, sub.phylo, tip.labels=T)
      subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
      
      #sub.populations is the subset of populations specific to a particular clade and timeslice
      sub.populations = subset(subset.populations, time.of.origin < t & time.of.extinction > t)
      
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
        
      output = rbind(output, cbind(sim=sim,clade.id = c, time = t, corr.results, gamma.stat = Gamma.stat,
                                     clade.richness = length(unique(sub.populations$spp.name))))
      print(c(c,t,date(),length(sub.clade.phylo$tip.label),extinct.time.most.recent));
        
    } # end sub clade for loop
    
  }; # end timeslice loop
  
  #write all of this output to files
  write.csv(output,paste(analysis_dir,"/SENC_Stats_sim",sim,".csv",sep=""),quote=F,row.names=F);
  analysis.end = date();
  #FIXME: store these warnings to a file, along with sim.id?
  #print(c(warnings(),sim.start,sim.end,analysis.end));
  
  
  ####################################################
  # Simulation summary plots
  ####################################################
  lat.grad.time.plot(sim.results, numslices = 10, output.dir = analysis_dir)
  clade.origin.corr.plot(output, params.out)
  clade.exmpl.figs(sim.results, output, clade.slices=6, seed=0, output.dir = analysis_dir)
  
  
} # end sim loop

