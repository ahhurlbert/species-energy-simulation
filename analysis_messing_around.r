setwd('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses')

####################
# Note: this code needs to be modified to take into account James' revised method of calculating clade origin times

require(ape)
require(geiger)
source('reg_calc_and_analysis_20121104.r') #added "integer" as acceptable class for max.time

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

par(mfrow=c(1,1))
plot(c(1,10),c(0,max(rich.time2$spp.rich)),type="n",xlab="Latitude",ylab="Species richness")
sapply(slices,function(x) points(11 - rich.time2$region[rich.time2$time==x], rich.time2$spp.rich[rich.time2$time==x],type='l',col=cols[which(slices==x)]))
legend('topright',legend=slices,lty='solid',col=cols[2:(numslices+1)])

#Node origin times based on phylogeny at end of simulation
all.dist <- dist.nodes(phylo.out)
root.dist <- all.dist[length(phylo.out$tip.label)+1, ]
rootdist = data.frame(clade.id = names(root.dist),RD = root.dist)

stats2 = merge(stats,rootdist,by='clade.id',all.x=T)

#Plotting
pdf(paste('corrs_vs_cladeage_sim',stats$sim[1],'.pdf',sep=''),height=6,width=9)
par(mfrow=c(2,3),oma=c(5,1,0,0),mar=c(2,5,2,2))
#for (i in c(611,764)) {
  t = max(stats2$time)                 #comment this out if looping over multiple time slices
  stats2x = subset(stats2, time==t)
  attach(stats2x)
  spline.df = 4
  plot(RD, r.env.rich, xlab="",ylab="Environment-Richness correlation")
  points(smooth.spline(RD[!is.na(r.env.rich)],r.env.rich[!is.na(r.env.rich)],df=spline.df),type='l',col='red')
  rect(-50,.5,t,1.1,col=rgb(.1,.1,.1,.1),lty=0)
  plot(RD, r.MRD.rich, xlab="",ylab="Mean Root Distance-Richness correlation")
  points(smooth.spline(RD[!is.na(r.MRD.rich)],r.MRD.rich[!is.na(r.MRD.rich)],df=spline.df),type='l',col='red')
  rect(-50,-1.1,t,-.5,col=rgb(.1,.1,.1,.1),lty=0)
  plot(RD, r.PSV.rich, xlab = "", ylab="PSV-Richness correlation")
  points(smooth.spline(RD[!is.na(r.PSV.rich)],r.PSV.rich[!is.na(r.PSV.rich)],df=spline.df),type='l',col='red')
  rect(-50,.5,t,1.1,col=rgb(.1,.1,.1,.1),lty=0)
  plot(RD, r.time.rich, xlab = "",ylab="Time in region-Richness correlation")
  points(smooth.spline(RD[!is.na(r.time.rich)],r.time.rich[!is.na(r.time.rich)],df=spline.df),type='l',col='red')
  rect(-50,.5,t,1.1,col=rgb(.1,.1,.1,.1),lty=0)
  plot(RD, r.env.MRD, xlab="", ylab="Environment-MRD correlation")
  points(smooth.spline(RD[!is.na(r.env.MRD)],r.env.MRD[!is.na(r.env.MRD)],df=spline.df),type='l',col='red')
  rect(-50,-1.1,t,-.5,col=rgb(.1,.1,.1,.1),lty=0)
  plot(RD, r.env.PSV, xlab="", ylab="Environment-PSV correlation")
  points(smooth.spline(RD[!is.na(r.env.PSV)],r.env.PSV[!is.na(r.env.PSV)],df=spline.df),type='l',col='red')
  rect(-50,.5,t,1.1,col=rgb(.1,.1,.1,.1),lty=0)
  mtext("Clade origin time",1,outer=T,line=2)
  detach(stats2x)
#}
dev.off()

#Looking at detailed clade level patterns for example clades of different age
clade.times = c(0,200,300,400,500,600)
#In line below, replace stats2$time with all.populations$time.of.origin once that gets fixed
clades = sapply(clade.times, function(x) {set.seed(1); sample(rootdist$clade.id[rootdist$RD> (x-10) & rootdist$RD < (x+10)],1)})
clds = as.numeric(as.character(clades))

pdf('clade_example_details.pdf',width=9,height=6)
par(mfrow=c(2,3),oma = c(5,1,4,0),mar=c(4,4,1,1))
for (i in clds) {
  cl.members = clade.members(i,phylo.out)
  cl.subset = subset(all.pops, spp.name %in% cl.members)
  cl.analysis = regional.calc(cl.subset[cl.subset<11,c('region','spp.name','time.of.origin')], phylo.out, max.time)
  attach(cl.analysis)
  plot(richness~region, xlab="Environment",ylab="Species Richness", cex=3, pch=16, col='darkblue', xlim=c(1,10))
  legend("top",legend=paste('r =',round(cor(richness,region,use="complete.obs"),2)), bty='n')
  plot(richness~MRD, xlab = "MRD", ylab="Species Richness", cex=3, pch=15, col = 'skyblue')
  legend("top",legend=paste('r =',round(cor(richness,MRD,use="complete.obs"),2)), bty='n')
  plot(richness~PSV, xlab = "PSV", ylab="Species Richness", cex=3, pch=17, col = 'darkgreen')
  legend("top",legend=paste('r =',round(cor(richness,PSV,use="complete.obs"),2)), bty='n')
  plot(richness~time.in.region, xlab = "Time in Region", ylab="Species Richness", cex=3, pch=17, col = 'yellow3')
  legend("top",legend=paste('r =',round(cor(richness,time.in.region,use="complete.obs"),2)), bty='n')
  plot(MRD~region, ylab = "MRD", xlab="Environment", cex=3, pch=18, col = 'orange', xlim=c(1,10))
  legend("top",legend=paste('r =',round(cor(region,MRD,use="complete.obs"),2)), bty='n')
  plot(PSV~region, ylab = "PSV", xlab="Environment", cex=3, pch=15, col = 'darkred', xlim=c(1,10))
  legend("top",legend=paste('r =',round(cor(region,PSV,use="complete.obs"),2)), bty='n')
  mtext(paste("CladeID =",i,"; Clade time of origin =", rootdist[rootdist$clade.id==i,'RD'],"; Clade richness =",length(cl.members)),3,outer=T)
  detach(cl.analysis)
}
dev.off()
