# Creating a supplemental figure depicting the phylogeny for one example simulation
# (simuation ID # 3325)

# First, read in simulation output
sim = 3325
stats = read.csv(paste(sim_dir,'/SENC_stats_sim',sim,'.csv',sep=''), header=T)
sim.results = output.unzip(sim_dir,sim)
all.populations = sim.results$all.populations
time.richness = sim.results$time.richness
phylo.out = sim.results$phylo.out
params.out = sim.results$params.out

# Evaluate phylogeny at equilibrium (t = 30000) and ignoring extinct species
t = 30000 
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

# Identify two example clades for highlighting
c2486 = clade.members(2486, sub.phylo, tip.labels=T)
c2486br = which.edge(sub.phylo, c2486)

c3157 = clade.members(3157, sub.phylo, tip.labels=T)
c3157br = which.edge(sub.phylo, c3157)

cols = rep('gray50', nrow(sub.phylo$edge))
cols[c2486br] = 'yellow'
cols[c3157br] = 'orange'
pdf('phylogeny_sim3325.pdf', height = 12, width = 12)
plot(sub.phylo, type='fan', show.tip.label=F, edge.col = cols)
dev.off()

