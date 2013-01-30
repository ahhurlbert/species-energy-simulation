# Sim 2635: trop origin, w=3, sigma=3, alpha=1e-04, specn = 1e-06, K gradient
library(caper)
library(ape)



sim_dir = 'c:/sencoutput/senc.out.130115'
sim = 2635

stats.output = read.csv(paste(sim_dir,'/SENC_Stats_sim',sim,'.csv',sep=''), header=T)
sim.out = output.unzip(sim_dir, sim)
phylo.out = sim.out$phylo.out
all.populations = sim.out$all.populations


#identify a clade with a strong negative richness-environment correlation that arises on the early side
plot(log10(stats.output$clade.origin.time), stats.output$r.env.rich,type="n")
text(log10(stats.output$clade.origin.time), stats.output$r.env.rich, stats.output$clade.id)

# Explore clade 3011
stats.output[stats.output$clade.id==3011,]

#From analysis_workflow_local.r
extant.ornot = aggregate(all.populations$extant,by=list(all.populations$spp.name),sum)
extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
t=30000
sub.species = as.character(unique(subset(all.populations,time.of.sp.origin <= t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);
sub.species2 = sub.species[!sub.species %in% extinct.species]
tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
sub.phylo = drop.tip(phylo.out,tips.to.drop);
mem.3011 = clade.members(3011,sub.phylo)
spp.name.3011=sub.phylo$tip.label[mem.3011]

# Find the node id in phylo.out that is equivalent to node 3011 in sub.phylo
phylo.3011 = drop.tip(phylo.out, as.character(phylo.out$tip.label[!phylo.out$tip.label %in% spp.name.3011]))
dist.3011 = cophenetic.phylo(phylo.3011) #distances between species
max.dist = max(dist.3011)
# Multiple distances will be equal to the max, but peruse the matrix to find two species that are this dist
# apart, e.g., 4885 and 3217

# calc most recent common ancestor
mrca.all = mrca(phylo.out) #takes awhile (e.g. >5 min, <1 hr)
mrca.3011=mrca.all["4885","3217"]
# we find that node # 11987 is the common ancestor

#verify by plotting
pdf('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/testphylo.pdf',height=72,width=72)
plot(phylo.out,type="fan")
tiplabels(spp.name.3011,which(phylo.out$tip.label %in% spp.name.3011))
nodelabels(11987,11987)
dev.off()

all.clade.members.3011 = clade.members(11987,phylo.out)
all.clade.member.names.3011 = phylo.out$tip.label[all.clade.members.3011]

clade.3011.pops = subset(all.populations, spp.name %in% all.clade.member.names.3011)