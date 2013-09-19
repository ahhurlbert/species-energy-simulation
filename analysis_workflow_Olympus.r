#!/usr/bin/env Rscript

sim = commandArgs();
sim = as.numeric(sim[length(sim)]);

# Choose number of time slices per simulation to analyze
num.of.time.slices = -999; # use -999 if you want to define specific time slices
which.time.slices = -999;
time.sequence = seq(2,100,by = 2);

# choose root only or all clades
root.only = 1 # 0 means all clades, 1 means just the root

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
  Rlib.location = "/pic/people/steg815/Rlibs"}

print(Rlib.location);
#print(length(list.files(Rlib.location)));

# Simulation workflow

#(2) load simulation and analysis functions
library(mnormt,lib.loc="/pic/people/steg815/Rlibs");
library(rgl,lib.loc="/pic/people/steg815/Rlibs");
library(ape,lib.loc="/pic/people/steg815/Rlibs");
library(permute,lib.loc="/pic/people/steg815/Rlibs");
library(nlme,lib.loc="/pic/people/steg815/Rlibs");
library(vegan,lib.loc="/pic/people/steg815/Rlibs");
library(picante,lib.loc="/pic/people/steg815/Rlibs");
library(mvtnorm,lib.loc="/pic/people/steg815/Rlibs");
library(caper,lib.loc="/pic/people/steg815/Rlibs");
library(paleotree,lib.loc="/pic/people/steg815/Rlibs");
library(plyr,lib.loc="/pic/people/steg815/Rlibs");
library(phytools, lib.loc="/pic/people/steg815/Rlibs");

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
sim.matrix$n.regions = NA
sim.matrix$extant.S = NA
sim.matrix$extinct.S = NA
sim.matrix$skipped.clades = NA
sim.matrix$skipped.times = NA
sim.matrix$BK.reg = NA
sim.matrix$BK.env = NA

trop.orig.extreme = 3;
temp.orig.extreme = 8;
pre.equil.time = 5459;

  output = numeric();
  
  # (5) read in simulation results for specified simulation from the output zip file
  sim.results = output.unzip(getwd(),sim)
  
  if ( !is.null(sim.results) ) {
    all.populations = sim.results$all.populations
    time.richness = sim.results$time.richness
    phylo.out = sim.results$phylo.out
    params.out = sim.results$params.out
    
    if (sim.matrix$reg.of.origin[sim.matrix$sim.id == sim] == 'tropical' ) { extreme.bin = trop.orig.extreme };
    if (sim.matrix$reg.of.origin[sim.matrix$sim.id == sim] == 'temperate' ) { extreme.bin = temp.orig.extreme };
    
    #if (time.richness$spp.rich[time.richness$time == pre.equil.time & time.richness$region == extreme.bin] > 5) {
    
    max.time.actual = max(time.richness$time);
    # If just a single timeslice, then use the end of the simulation or a designated time, otherwise space them equally
    if (num.of.time.slices == 1) { timeslices = pre.equil.time};
    if (which.time.slices != -999  ) { timeslices = which.time.slices };
    if (num.of.time.slices > 1) {timeslices = as.integer(round(seq(max(time.richness$time)/num.of.time.slices,max(time.richness$time),length=num.of.time.slices),digits=0))};
    if (time.sequence[1] != -999) {timeslices = subset(time.sequence,time.sequence <= max(time.richness$time))};
    
    skipped.clades = 0
    skipped.times = ""
    for (t in timeslices) {
      # vector of species in existence at time t
      sub.species = as.character(unique(subset(all.populations,time.of.sp.origin <= t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);
      
      # Some species may be extant globally (extant==1) but in our boundary regions (0,11) only;
      # we need to eliminate species that are not extant within regions 1-10 (which is all that is
      # reflected in the all.populations dataframe)
      time.slice.populations = all.populations;
      time.slice.populations$extant = 0;
      time.slice.populations$extant[time.slice.populations$time.of.origin <= t & time.slice.populations$time.of.extinction > t] = 1
      extant.ornot = aggregate(time.slice.populations$extant,by=list(time.slice.populations$spp.name),sum)
      extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
      
      # FIXME:
      # Add more explanatory comments justifying why we don't need to consider species that existed
      # at time t but went extinct before the present.
      # In some cases (e.g. sim 1 or 2, t=6000), tips.to.drop includes all tips and so sub.phylo is empty.
      # Does it make sense for this to ever happen? If not, fix it.
      # If so, need to provide an if-else error catch both in the creation of sub.phylo,
      # and of sub.clade.phylo inside the clade loop. (Sim 3, t = 156 bonks at that point)
      # NOTE: code runs for sim==5 currently as a test case
      sub.species2 = sub.species[!sub.species %in% extinct.species]
      tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
      
      # check to see if there are at least min.num.spp species for continuing with the analysis; if not store the skipped timeslice
      if ( (length(phylo.out$tip.label) - length(tips.to.drop)) < min.num.spp) {
        skipped.times = paste(skipped.times, t) # keep track of the timeslices that were skipped in a text string
      } else {
        
        sub.phylo = drop.tip(phylo.out,tips.to.drop);
        temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1]); temp.root.time;
        most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1])))]; most.recent.spp;
        extinct.time.most.recent = unique(all.populations$time.of.sp.extinction[all.populations$spp.name==most.recent.spp]); extinct.time.most.recent;
        sub.phylo$root.time = temp.root.time + max(c(0,max.time.actual - extinct.time.most.recent)); sub.phylo$root.time;
        sub.phylo = collapse.singles(timeSliceTree(sub.phylo,sliceTime=(max.time.actual - t),plot=F,drop.extinct = T));
        num.of.spp = length(sub.phylo$tip.label);
        
        if (root.only == 1) { sub.clade.loop.end = (num.of.spp+1) }
        if (root.only == 0) { sub.clade.loop.end = max(sub.phylo$edge) }
        for (c in (num.of.spp+1):sub.clade.loop.end) {
          
          #pull out list of species names belonging to each subclade
          sub.clade = clade.members(c, sub.phylo, tip.labels=T)
          subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
          
          #sub.populations is the subset of populations specific to a particular clade and timeslice
          sub.populations = subset(subset.populations, time.of.origin <= t & time.of.extinction > t)
          
          #sub.clade.phylo is a specific simulation clade pulled from the phylogeny that was sliced at timeslice t
          tips.to.drop2 = as.character(sub.phylo$tip.label[which(is.element(sub.phylo$tip.label,as.character(sub.populations$spp.name))==F)]);
          
          # check to see if there are at least min.num.spp species for continuing with the analysis; if not increment skipped.clades
          if((length(sub.phylo$tip.label) - length(tips.to.drop2)) < min.num.spp) {
            skipped.clades = skipped.clades + 1
          } else {
            
            sub.clade.phylo = drop.tip(sub.phylo,tips.to.drop2);
            sub.clade.phylo$root.time = max(dist.nodes(sub.clade.phylo)[1:Ntip(sub.clade.phylo),Ntip(sub.clade.phylo) + 1]); sub.clade.phylo$root.time;
            sub.clade.phylo$origin.time = t - sub.clade.phylo$root.time; sub.clade.phylo$origin.time;
            
            if (identical(sort(as.integer(unique(sub.populations$spp.name))) , sort(as.integer(sub.clade.phylo$tip.label)))==F ) {print(c(c,t,'Error: trimmed phylogeny does not contain the correct species')); break} else{};
            
            reg.summary = regional.calc(sub.populations[,c('region','spp.name','time.of.origin','reg.env','extant')], sub.clade.phylo, as.integer(t));
            
            #Note that extinction calculation must be done on subset.populations, not sub.populations
            extinction = extinct.calc(subset.populations, timeslice=t)
            reg.summary2 = merge(reg.summary,extinction[,c('region','extinction.rate')],by='region')
            
            corr.results = cbind(xregion.analysis(reg.summary2),length(reg.summary$region[reg.summary$richness > 1]))
            colnames(corr.results)[ncol(corr.results)] = 'n.div.regions'
            
            #Pybus & Harvey (2000)'s gamma statistic
            Gamma.stat = gammaStat(sub.clade.phylo)
            
            #Calculate Blomberg's K for two traits: environmental optimum, and mean region of occurrence
            spp.traits = aggregate(sub.populations$region, by = list(sub.populations$spp.name, sub.populations$env.opt),
                                   function(x) mean(x, na.rm=T))
            names(spp.traits) = c('spp.name','env.opt','region')
            
            spp.env = spp.traits$env.opt
            names(spp.env) = spp.traits$spp.name
            BK.env = phylosig(sub.clade.phylo, spp.env[sub.clade.phylo$tip.label], method="K")
            
            spp.reg = spp.traits$region
            names(spp.reg) = spp.traits$spp.name
            BK.reg = phylosig(sub.clade.phylo, spp.reg[sub.clade.phylo$tip.label], method="K")
            
            output = rbind(output, cbind(sim=sim,clade.id = c, time = t, corr.results, gamma.stat = Gamma.stat,
                                         clade.richness = length(unique(sub.populations$spp.name)), 
                                         BK.env = BK.env , BK.reg = BK.reg))
            print(paste(sim,sub.clade.loop.end,c,t,date(),length(sub.clade.phylo$tip.label),sep="   "));
          } # end third else
        } # end sub clade for loop
      } # end second else
    }; # end timeslice loop
    
    #write all of this output to files
    if (num.of.time.slices == 1) {write.csv(output,paste("NEW_Stats_sim",sim,"_time",t,".csv",sep=""),quote=F,row.names=F)};
    if (num.of.time.slices > 1 & root.only == 0) {write.csv(output,paste("NEW_Stats_sim",sim,"_mult_times_all_clades.csv",sep=""),quote=F,row.names=F)};
    if (num.of.time.slices > 1 & root.only == 1) {write.csv(output,paste("NEW_Stats_sim",sim,"_mult_times_root_only.csv",sep=""),quote=F,row.names=F)};
    if (which.time.slices != -999 & root.only == 1) {write.csv(output,paste("NEW_Stats_sim",sim,"_specific_times_root_only.csv",sep=""),quote=F,row.names=F)};
    if (which.time.slices != -999 & root.only == 0) {write.csv(output,paste("NEW_Stats_sim",sim,"_specific_times_all_clades.csv",sep=""),quote=F,row.names=F)};
    if (time.sequence != -999 & root.only == 1) {write.csv(output,paste("NEW_Stats_sim",sim,"_time_seq_root_only.csv",sep=""),quote=F,row.names=F)};
    if (time.sequence != -999 & root.only == 0) {write.csv(output,paste("NEW_Stats_sim",sim,"_time_seq_all_clades.csv",sep=""),quote=F,row.names=F)};
    
    analysis.end = date();
    #FIXME: store these warnings to a file, along with sim.id? Or is this being done in the shell?
    #print(c(warnings(),sim.start,sim.end,analysis.end));
    
    
    ####################################################
    # Simulation summary plots
    ####################################################
    lat.grad.time.plot(sim.results, numslices = 10, output.dir = getwd())
    
    # clade.origin.corr.plot only if there are some output rows 
    if(!is.null(nrow(output))) {
      clade.origin.corr.plot(output, params.out, output.dir = getwd())
      
      # Number of rows of output with at least 1 correlation (there are 4 non-correlation cols in corr.results)
      sim.matrix[sim.matrix$sim.id==sim,'output.rows'] = sum(apply(output,1,function(x) sum(is.na(x)) < (ncol(corr.results)-4)))
    } else {
      sim.matrix[sim.matrix$sim.id==sim,'output.rows'] = 0
    }
    
    # There are currently some bugs in clade.exmpl.figs.
    # clade.exmpl.figs(sim.results, output, clade.slices=6, seed=0, output.dir = analysis_dir)
    
    # Add overall summary info
    sim.matrix[sim.matrix$sim.id==sim,'n.regions'] = length(unique(all.populations$region))
    sim.matrix[sim.matrix$sim.id==sim,'extant.S'] = nrow(extant.ornot[extant.ornot$x>0,])
    sim.matrix[sim.matrix$sim.id==sim,'extinct.S'] = length(extinct.species)
    sim.matrix[sim.matrix$sim.id==sim,'skipped.clades'] = skipped.clades # number of clades skipped over for analysis, summed over timeslices
    sim.matrix[sim.matrix$sim.id==sim,'skipped.times'] = skipped.times # number of time slices skipped over for analysis
    sim.matrix[sim.matrix$sim.id==sim,'BK.reg'] = BK.reg # blomberg's K based on region
    sim.matrix[sim.matrix$sim.id==sim,'BK.env'] = BK.env # blomberg's K based on environment
    
    write.csv(sim.matrix[sim.matrix$sim.id==sim,],paste("sim.matrix.output.",sim,"_time",t,".csv",sep=""),quote=F,row.names=F);
    sim.matrix[sim.matrix$sim.id==sim,]
    #} else {print(sim)} # end second if (richness check for defined, extreme bin)
  } # end first if (file check)
