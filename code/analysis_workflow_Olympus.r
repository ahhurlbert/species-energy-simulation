#!/usr/bin/env Rscript

if (local == 0) {
  sim = commandArgs();
  sim = as.numeric(sim[length(sim)]);
}
  
# Choose number of time slices per simulation to analyze
num.of.time.slices = 1; # use -999 if you want to define specific time slices

# which.time.slices is (apparently) only for specifying particular, unevenly spaced time slices;
# if not being used it should be set to -999
which.time.slices = -999;

# time.sequence is (apparently) for when the time slices occur for a regular interval; set to -999 if not being used
# Note that due to the slow calculation of tree imbalance (beta) for large trees, it may be best to specify only ~20 time slices
time.sequence = -999; # to turn off 'time.sequence'
#time.sequence = seq(2,300,by=2); # for time scenario sims
#time.sequence = seq(1000,100000,length=100); # for energy gradient sims

# choose root only or all clades
root.only = 0 # 0 means all clades, 1 means just the root

# Set minimum number of species in a clade needed to proceed with analysis
min.num.spp = 8;

# Specify number of clusters available for parallel proccessing
num.clusters = 2

# Simulation workflow

#(2) load simulation and analysis functions
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

package.vector = c('ape','permute','nlme','vegan','picante','mvtnorm','caper','paleotree','plyr','phytools','apTreeshape');

source('code/reg_calc_and_analysis.r');
source('code/make.phylo.jimmy.fun.r');
source('code/lat.grad.time.plot.r');
source('code/clade.origin.corr.plot.r');
source('code/clade.exmpl.figs.r');
source('code/extinct.calc.r');
source('code/unzipping_files.r');
source('code/maxlik.betasplit.AH.r');

if (local == 1) {
  cl = makeCluster(num.clusters);
  registerDoParallel(cl);
}

#(3) read in master simulation matrix with chosen parameter combinations;
# then add fields for storing output summary
sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv",header=T)
sim.matrix$n.regions = NA
sim.matrix$extant.S = NA
sim.matrix$extinct.S = NA
sim.matrix$skipped.clades = NA
sim.matrix$skipped.times = NA

  rm(list=c('all.populations', 'time.richness', 'phylo.out', 'params.osput', 'output', 'sim.results'))
  output = numeric();
  
  # (5) read in simulation results for specified simulation from the output zip file, or the raw output files
  sim.results = output.unzip('raw_sim_output', sim)
  
  if ( !is.null(sim.results) ) {
    all.populations = sim.results$all.populations
    time.richness = sim.results$time.richness
    phylo.out = sim.results$phylo.out
    params.out = sim.results$params.out
    
    max.time.actual = max(time.richness$time);
    # If just a single timeslice, then use the end of the simulation or a designated time, otherwise space them equally (which.time.slices == -999)
    # or use specified vector in which.time.slices
    if (num.of.time.slices == 1) { 
      timeslices = max.time.actual 
    } else {
      if (which.time.slices != -999 & num.of.time.slices == - 999) { timeslices = which.time.slices };
      if (which.time.slices == -999 & num.of.time.slices > 1) {
        timeslices = as.integer(round(seq(max(time.richness$time) / num.of.time.slices,
                                          max(time.richness$time), length = num.of.time.slices), digits = 0))
      }
      if (time.sequence[1] != -999) { timeslices = subset(time.sequence, time.sequence <= max.time.actual) }
    }
    
    # Conduct analysis for each time slice specified above
    skipped.clades = 0
    skipped.times = ""
    for (t in timeslices) {
      # vector of species in existence at time t
      sub.species = as.character(unique(subset(all.populations, time.of.origin <= t & time.of.extinction > t, select = 'spp.name'))[,1]);
      
      # Some species may be extant globally (extant==1) but in our boundary regions (0 or 11) only;
      # we need to eliminate species that are not extant within regions 1-10 (which is all that is
      # reflected in the all.populations dataframe)
      time.slice.populations = all.populations;
      time.slice.populations$extant = 0;
      time.slice.populations$extant[time.slice.populations$time.of.origin <= t & time.slice.populations$time.of.extinction > t] = 1
      extant.ornot = aggregate(time.slice.populations$extant, by = list(time.slice.populations$spp.name), sum)
      extinct.species = as.character(extant.ornot[extant.ornot$x == 0, 'Group.1'])
      
      sub.species2 = sub.species[!sub.species %in% extinct.species] # 'true' list of extant species
      tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]) # species to drop
      
      # check to see if there are at least min.num.spp species for continuing with the analysis; if not store the skipped timeslice
      if ( (length(phylo.out$tip.label) - length(tips.to.drop)) < min.num.spp) {
        skipped.times = paste(skipped.times, t) # keep track of the timeslices that were skipped in a text string
      } else {
        
        # Drop extinct species out of the phylogeny and tidy up
        sub.phylo = drop.tip(phylo.out, tips.to.drop);
        temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo), Ntip(sub.phylo) + 1])
        most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo), Ntip(sub.phylo) + 1])))]
        extinct.time.most.recent = unique(all.populations$time.of.sp.extinction[all.populations$spp.name == most.recent.spp])
        sub.phylo$root.time = temp.root.time + max(c(0, max.time.actual - extinct.time.most.recent))
        sub.phylo = collapse.singles(timeSliceTree(sub.phylo, sliceTime = (max.time.actual - t), plot = F, drop.extinct = T));
        num.of.spp = length(sub.phylo$tip.label);
        
        
        # The 'c' loop repeats the basic analysis for every node in the tree (if root.only = 0)
        # or just for the overall root clade (if root.only = 1)
        if (root.only == 1) { sub.clade.loop.end = (num.of.spp+1) }
        if (root.only == 0) { sub.clade.loop.end = max(sub.phylo$edge) }
        for (c in (num.of.spp+1):sub.clade.loop.end) {
          
          #pull out list of species names belonging to each subclade
          sub.clade = clade.members(c, sub.phylo, tip.labels=T)
          subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
          
          #sub.populations is the subset of populations specific to a particular clade and timeslice
          sub.populations = subset(subset.populations, time.of.origin <= t & time.of.extinction > t)
          
          #sub.clade.phylo is a specific simulation clade pulled from the phylogeny that was sliced at timeslice t
          tips.to.drop2 = as.character(sub.phylo$tip.label[which(is.element(sub.phylo$tip.label, as.character(sub.populations$spp.name)) == F)]);
          
          # check to see if there are at least min.num.spp species for continuing with the analysis; if not increment skipped.clades
          if((length(sub.phylo$tip.label) - length(tips.to.drop2)) < min.num.spp) {
            skipped.clades = skipped.clades + 1
          } else {
            
            sub.clade.phylo = drop.tip(sub.phylo, tips.to.drop2);
            sub.clade.phylo$root.time = max(dist.nodes(sub.clade.phylo)[1:Ntip(sub.clade.phylo), Ntip(sub.clade.phylo) + 1])
            sub.clade.phylo$origin.time = t - sub.clade.phylo$root.time
            
            if (identical(sort(as.integer(unique(sub.populations$spp.name))) , sort(as.integer(sub.clade.phylo$tip.label))) == F ) {
              print(c(c, t, 'Error: trimmed phylogeny does not contain the correct species'))
              break
            } 
            
            # Calculate summary statistics at the regional level
            reg.summary = regional.calc(sub.populations[,c('region','spp.name','time.of.origin','reg.env','extant')], 
                                        sub.clade.phylo, 
                                        as.integer(t));
            
            #Note that extinction calculation must be done on subset.populations, not sub.populations
            extinction = extinct.calc(subset.populations, timeslice = t)
            reg.summary2 = merge(reg.summary, extinction[, c('region', 'extinction.rate')], by = 'region')
            
            MRD.range = max(reg.summary$MRD,na.rm = T) - min(reg.summary$MRD,na.rm = T)
            MRD.mean = mean(reg.summary$MRD,na.rm = T)
            MRD.var = var(reg.summary$MRD,na.rm = T)
            MRD.rich.slope = lm(reg.summary$MRD ~ reg.summary$richness)$coefficients[2]
            MRD.env.slope = lm(reg.summary$MRD ~ reg.summary$reg.env)$coefficients[2]
            PSV.range = max(reg.summary$PSV,na.rm = T) - min(reg.summary$PSV,na.rm = T)
            PSV.mean = mean(reg.summary$PSV,na.rm = T)
            PSV.var = var(reg.summary$PSV,na.rm = T)
            PSV.rich.slope = lm(reg.summary$PSV ~ reg.summary$richness)$coefficients[2]
            PSV.env.slope = lm(reg.summary$PSV ~ reg.summary$reg.env)$coefficients[2]
            n.div.regions = length(reg.summary$region[reg.summary$richness > 1])
            
            corr.results = cbind(xregion.analysis(reg.summary2),
                                 MRD.range,
                                 MRD.mean,
                                 MRD.var,
                                 MRD.rich.slope,
                                 MRD.env.slope,
                                 PSV.range,
                                 PSV.mean,
                                 PSV.var,
                                 PSV.rich.slope,
                                 PSV.env.slope,
                                 n.div.regions)
            
            #Pybus & Harvey (2000)'s gamma statistic
            Gamma.stat = gammaStat(sub.clade.phylo)
            
            #Calculate Blum & Francois (2006)'s Beta metric of tree imbalance using apTreeshape package
            # --seems to bonk on very large phylogenies, so only try calculating for fewer than 6000 species
            if(length(sub.phylo$tip.label) < 100000) {
              tree.beta.out = maxlik.betasplit.AH(sub.clade.phylo)
              tree.beta = tree.beta.out$max_lik
            } else {
              tree.beta = NA
            }
            
            output = rbind(output, cbind(sim = sim, clade.id = c, time = t, corr.results, gamma.stat = Gamma.stat,
                                         clade.richness = length(unique(sub.populations$spp.name)), 
                                         tree.beta = tree.beta))
            print(paste(sim, sub.clade.loop.end, c, t, date(), length(sub.clade.phylo$tip.label), sep="   "));
            flush.console();
          } # end third else
        } # end sub clade for loop
      } # end second else
    }; # end timeslice loop

    # this print statement should show up in the output files on the cluster.
    print(warnings())

    #write all of this output to files
    if (root.only == 0) {
      write.csv(output, paste("Stats_sim", sim, "_all_subclades.csv", sep = ""), quote = F, row.names = F)
    }
    if (root.only == 1) {
      write.csv(output, paste("Stats_sim", sim, "_rootclade_only.csv", sep = ""), quote = F, row.names = F)
    }
    
    analysis.end = date();
    
    # Add overall summary info
    sim.matrix[sim.matrix$sim.id == sim, 'n.regions'] = length(unique(all.populations$region))
    sim.matrix[sim.matrix$sim.id == sim, 'extant.S'] = nrow(extant.ornot[extant.ornot$x>0,])
    sim.matrix[sim.matrix$sim.id == sim, 'extinct.S'] = length(extinct.species)
    sim.matrix[sim.matrix$sim.id == sim, 'skipped.clades'] = skipped.clades # number of clades skipped over for analysis, summed over timeslices
    sim.matrix[sim.matrix$sim.id == sim, 'skipped.times'] = skipped.times # number of time slices skipped over for analysis
    
    write.csv(sim.matrix[sim.matrix$sim.id == sim,], paste("analysis_output/summary_output_sim", sim, ".csv", sep = ""), quote = F, row.names = F);
  } # end first if (file check)

