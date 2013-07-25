source('unzipping_files.r')

library(picante);
library(mvtnorm);
library(caper);
library(paleotree);
library(plyr);
library(phytools);


# This function compiles the following info across the latitudinal gradient for a particular clade (cladeid):
# --overal equilibrial richness
# --clade richness
# --time of colonization of each region
# --number of species in each region at the respective time of colonization
# When output=T, a dataframe with this info is produced; otherwise, only a graph of these values
reg.opp.grad = function(cladeid, sim, sim_dir = 'C:/SENCoutput', t = 30000, output = T, log=T, plot=F) {
  stats = read.csv(paste(sim_dir,'/SENC_stats_sim',sim,'.csv',sep=''), header=T)
  sim.results = output.unzip(sim_dir,sim)
  all.populations = sim.results$all.populations
  time.richness = sim.results$time.richness
  phylo.out = sim.results$phylo.out
  params.out = sim.results$params.out
  
  #############################################
  # Code chunk from 'analysis_workflow_local.r'
  #############################################
  t = 30000 #richness evaluated at the end of the simulation
  max.time.actual = max(time.richness$time)
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
  sub.species2 = sub.species[!sub.species %in% extinct.species]
  tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
  sub.phylo = drop.tip(phylo.out,tips.to.drop);
  temp.root.time = max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1]); temp.root.time;
  most.recent.spp = sub.phylo$tip.label[as.numeric(names(which.max(dist.nodes(sub.phylo)[1:Ntip(sub.phylo),Ntip(sub.phylo) + 1])))]; most.recent.spp;
  extinct.time.most.recent = unique(all.populations$time.of.sp.extinction[all.populations$spp.name==most.recent.spp]); extinct.time.most.recent;
  sub.phylo$root.time = temp.root.time + max(c(0,max.time.actual - extinct.time.most.recent)); sub.phylo$root.time;
  sub.phylo = collapse.singles(timeSliceTree(sub.phylo,sliceTime=(max.time.actual - t),plot=F,drop.extinct = T));
  num.of.spp = length(sub.phylo$tip.label);
  #############################################
  
  
  sub.clade = clade.members(cladeid, sub.phylo, tip.labels=T)
  subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
  sub.populations = subset(subset.populations, time.of.origin <= t & time.of.extinction > t)
  clade.reg.rich = data.frame(table(sub.populations$region))
  names(clade.reg.rich) = c('region','clade.rich')
  #Richness gradient at the end of the simulation (t=30000)
  total.latgrad = subset(time.richness, time==30000 & region > 0 & region < 11, select = c('region','spp.rich'))
  #Carrying capacity (individuals) across the gradient
  total.latgrad$K = seq(4000, 40000, length.out = 10)
  regions = unique(subset.populations$region)
  reg.origin.time = aggregate(subset.populations$time.of.origin, by = list(subset.populations$region), min)
  names(reg.origin.time) = c('region','time')
  reg.origin.time$time[reg.origin.time$time==0] <- 1
  reg.origin.time$rich.at.time = sapply(1:nrow(reg.origin.time), function(x) 
    time.richness$spp.rich[time.richness$region==reg.origin.time$region[x] & time.richness$time==reg.origin.time$time[x]])
  #adding 1 in the function above so that if no species are present when clade colonizes, richness is now 1, etc.
  out = merge(total.latgrad, reg.origin.time, by = 'region', all.x=T)
  out$opp.frac = 1 - out$rich.at.time/out$spp.rich   # Opportunity as measured by the fraction of equilibrial richness not yet present at time of colonization
  out$opp.diff = out$spp.rich - out$rich.at.time     # Opportunity as measured by the difference between equilibrial richness and actual richness at time of colonization
  out$opp.frac[out$opp.frac < 0] <- 0
  out$opp.diff[out$opp.diff < 0] <- 0
  out$opp.K = out$K/out$rich.at.time      # Opportunity as measured by potential mean population size at time of colonization
  out2 = merge(out, clade.reg.rich, by = 'region', all.x=T)
  reg.of.origin = sub.populations$region[sub.populations$time.of.sp.origin == min(sub.populations$time.of.sp.origin)]
  out2$reg.of.origin = reg.of.origin
  
  if(plot) {
    if(log) {
    plot(11-out2$region,log10(out2$spp.rich), type='l', xlim=c(1,10), ylim = c(0,3),
         xlab = 'Latitudinal bin', ylab = 'log10 Species richness', 
         main = paste('Clade',cladeid, ', richness =',stats$clade.richness[stats$clade.id==cladeid], 
                      ', time of origin =', min(reg.origin.time$time)))
    points(11-out2$region, log10(out2$rich.at.time), type='l', col='red')
    points(11-out2$region, log10(out2$clade.rich), type = 'l', col='blue')
    arrows(11-reg.of.origin, 0.25*max(log10(out2$spp.rich)), 11-reg.of.origin, 0)
    par(new=T)
    plot(11-out2$region, log10(out2$opp.diff+1), type = 'l', xlim = c(1,10), col='green', yaxt="n", xlab="", ylab="", ylim = c(0,3))
    }
    else {
      plot(11-out2$region,out2$spp.rich, type='l', xlim=c(1,10), ylim = c(0,500),
           xlab = 'Latitudinal bin', ylab = 'log10 Species richness', 
           main = paste('Clade',cladeid, ', richness =',stats$clade.richness[stats$clade.id==cladeid], 
                        ', time of origin =', min(reg.origin.time$time)))
      points(11-out2$region, out2$rich.at.time, type='l', col='red')
      points(11-out2$region, out2$clade.rich, type = 'l', col='blue')
      arrows(11-reg.of.origin, 0.25*max(out2$spp.rich), 11-reg.of.origin, 0)
      par(new=T)
      plot(11-out2$region, out2$opp.frac, type = 'l', xlim = c(1,10), col='green', yaxt="n", xlab="", ylab="", ylim=c(0,1))
    }
    axis(4)
    mtext("Proportion of equilibrial richness available", 4, line = 3)
  
    legend('topright',lty='solid',col = c('black','blue','red','green'), legend = c('Total richness','Clade richness',
                                                                                  'Richness at time\nof colonization',
                                                                                  'Mean population size at\ntime of colonization'))
  } #end if plot
  if(output) {
    return(out2)
  }
  
}

#pdf('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/richness_opportunity_gradients_RichFrac.pdf', height = 6, width =8)
#sapply(stats.3325$clade.id[stats.3325$n.regions > 3], function(x) reg.opp.grad(x, stats.3325, output=F, log=F))
#dev.off()


# From inspection of the plots, we settled on clade 3157 which exhibits a reverse latitudinal gradient,
# and clade 2486 which exhibits a fairly flat/humped shaped gradient



