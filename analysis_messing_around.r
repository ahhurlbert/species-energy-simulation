setwd('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses')

require(ape)
require(geiger)

phylo.out = read.tree('SENC_phylo_sim316.tre')
all.pops = read.csv('SENC_all.pops_sim316.csv',header=T)
rich.time = read.csv('SENC_time.rich_sim316.csv',header=T)
params = read.csv('SENC_params.out_sim316.csv',header=T)
stats = read.csv('SENC_stats_sim316.csv',header=T)

rich.time2 = rich.time[!rich.time$region %in% c(0,11),]

max.time = max(rich.time$time)
numslices = 10
slices = round(seq(max.time/numslices, max.time,by=max.time/numslices),0)
cols = rainbow(numslices+1)

plot(c(1,10),c(0,max(rich.time2$spp.rich)),type="n",xlab="Latitude",ylab="Species richness")
sapply(slices,function(x) points(11 - rich.time2$region[rich.time2$time==x], rich.time2$spp.rich[rich.time2$time==x],type='l',col=cols[which(slices==x)]))

all.dist <- dist.nodes(phylo.out)
root.dist <- all.dist[length(phylo.out$tip.label)+1, ]

#minor change here