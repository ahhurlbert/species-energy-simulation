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

#Node origin times based on phylogeny at end of simulation
all.dist <- dist.nodes(phylo.out)
root.dist <- all.dist[length(phylo.out$tip.label)+1, ]
rootdist = data.frame(clade.id = names(root.dist),RD = root.dist)

stats2 = merge(stats,rootdist,by='clade.id',all.x=T)

#Plotting
pdf()
par(mfrow=c(2,3),oma=c(5,1,0,0),mar=c(2,5,2,2))
attach(stats2)
plot(RD, r.env.rich, xlab="",ylab="Environment-Richness correlation")
points(smooth.spline(RD[!is.na(r.env.rich)],r.env.rich[!is.na(r.env.rich)],df=4),type='l',col='red')
plot(RD, r.MRD.rich, xlab="",ylab="Mean Root Distance-Richness correlation")
points(smooth.spline(RD[!is.na(r.MRD.rich)],r.MRD.rich[!is.na(r.MRD.rich)],df=4),type='l',col='red')
plot(RD, r.PSV.rich, xlab = "", ylab="PSV-Richness correlation")
points(smooth.spline(RD[!is.na(r.PSV.rich)],r.PSV.rich[!is.na(r.PSV.rich)],df=4),type='l',col='red')
plot(RD, r.time.rich, xlab = "",ylab="Time in region-Richness correlation")
points(smooth.spline(RD[!is.na(r.time.rich)],r.time.rich[!is.na(r.time.rich)],df=4),type='l',col='red')
plot(RD, r.env.MRD, xlab="", ylab="Environment-MRD correlation")
points(smooth.spline(RD[!is.na(r.env.MRD)],r.env.MRD[!is.na(r.env.MRD)],df=4),type='l',col='red')
plot(RD, r.env.PSV, xlab="", ylab="Environment-PSV correlation")
points(smooth.spline(RD[!is.na(r.env.PSV)],r.env.PSV[!is.na(r.env.PSV)],df=4),type='l',col='red')
mtext("Clade origin time",1,outer=T,line=2)
