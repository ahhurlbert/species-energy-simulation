# Sim 2635: trop origin, w=3, sigma=3, alpha=1e-04, specn = 1e-06, K gradient
user='Allen'
if(user=='Allen') {
  repo_dir = 'C:/Documents and Settings/Hurlbert/species-energy-simulation'
  sim_dir = 'c:/sencoutput/senc.out.130115'
  newplot_dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses'
  
}

library(caper)
library(ape)
source(paste(repo_dir,'/unzipping_files.r',sep=''))



#Note: This code currently seems to work for the sim/clade #s below, however,
#      for sim 3365, most clades have far fewer extant species than are
#      indicated in the Stats file for clade.richness.
#      Possibly these missing species are all in the edge regions (0,11)?

sim = 2635

stats.output = read.csv(paste(sim_dir,'/SENC_Stats_sim',sim,'.csv',sep=''), header=T)
sim.out = output.unzip(sim_dir, sim)
phylo.out = sim.out$phylo.out
all.populations = sim.out$all.populations


#identify a clade with a strong negative richness-environment correlation that arises on the early side
plot(log10(stats.output$clade.origin.time), stats.output$r.env.rich,type="n")
text(log10(stats.output$clade.origin.time), stats.output$r.env.rich, stats.output$clade.id)

# Explore clade 3011; note this is the node ID on sub.phylo, not phylo.out
clade = 3011
stats.output[stats.output$clade.id==clade,]

#From analysis_workflow_local.r
extant.ornot = aggregate(all.populations$extant,by=list(all.populations$spp.name),sum)
extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
t=30000
sub.species = as.character(unique(subset(all.populations,time.of.sp.origin <= t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);
sub.species2 = sub.species[!sub.species %in% extinct.species]
tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
sub.phylo = drop.tip(phylo.out,tips.to.drop);
mem.3011 = clade.members(clade,sub.phylo)
spp.name.3011=sub.phylo$tip.label[mem.3011]

# Find the node id in phylo.out that is equivalent to node 3011 in sub.phylo
phylo.3011 = drop.tip(phylo.out, as.character(phylo.out$tip.label[!phylo.out$tip.label %in% spp.name.3011]))
dist.3011 = cophenetic.phylo(phylo.3011) #distances between species
max.dist = max(dist.3011)
# Multiple distances will be equal to the max, but peruse the matrix to find two species that are this dist
# apart, e.g., 4885 and 3217
dist.3011

# calc most recent common ancestor
mrca.all = mrca(phylo.out) #takes awhile (e.g. >5 min, <1 hr)
mrca.3011=mrca.all["4885","3217"]
# we find that node # 11987 is the common ancestor
all.clade.members.3011 = clade.members(11987,phylo.out)
all.clade.member.names.3011 = phylo.out$tip.label[all.clade.members.3011]

#verify by plotting
pdf('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/testphylo.pdf',height=72,width=72)
plot(phylo.out,type="fan")
tiplabels(all.clade.member.names.3011,all.clade.members.3011,bg='green')
tiplabels(spp.name.3011,which(phylo.out$tip.label %in% spp.name.3011))
nodelabels(11987,11987)
dev.off()

#Subset all.populations to only include these clade members
clade.3011.pops = subset(all.populations, spp.name %in% all.clade.member.names.3011)

#timeslices = seq(600,5400,length.out = 9)     #set appropriately for the specific clade
time.window = timeslices[2] - timeslices[1]

reg.rich.thru.time = data.frame(time=NA, 
                                region=NA, 
                                total.rich=NA,       # richness of all species extant during the timeslice
                                clade.rich=NA,       # richness of clade members extant during the timeslice
                                total.ext.pops=NA,   # total number of populations going extinct in the window prior to timeslice
                                clade.ext.pops=NA,   # total number of pops within the clade going extinct prior to timeslice
                                total.ext.spp=NA,    # total number of species going globally extinct in the time window
                                clade.ext.spp=NA)    # total number of species within the clade going globally extinct in the time window

for (t in timeslices) {
  sub.pops = subset(clade.3011.pops, time.of.origin < t & time.of.extinction > t)
  cl.reg.rich = data.frame(table(sub.pops$region))

  cl.ext.pops = subset(clade.3011.pops, time.of.extinction <= t & time.of.extinction > (t - time.window))
  cl.reg.ext.pops = data.frame(table(cl.ext.pops$region))
  if (nrow(cl.reg.ext.pops)==0) { cl.reg.ext.pops = data.frame(Var1=1:10, Freq=rep(NA,10))}
  cl.ext.spp = subset(clade.3011.pops, time.of.sp.extinction <=t & time.of.sp.extinction > (t- time.window))
  cl.reg.ext.spp = data.frame(table(cl.ext.spp$region))
  if (nrow(cl.reg.ext.spp)==0) { cl.reg.ext.spp = data.frame(Var1=1:10, Freq=rep(NA,10))}
  
  all.pops = subset(all.populations, time.of.origin < t & time.of.extinction > t)
  all.reg.rich = data.frame(table(all.pops$region))
  
  all.ext.pops = subset(all.populations, time.of.extinction <= t & time.of.extinction > (t - time.window))
  all.reg.ext.pops = data.frame(table(all.ext.pops$region))
  if (nrow(all.reg.ext.pops)==0) { all.reg.ext.pops = data.frame(Var1=1:10, Freq=rep(NA,10))}
  all.ext.spp = subset(all.populations, time.of.sp.extinction <=t & time.of.sp.extinction > (t- time.window))
  all.reg.ext.spp = data.frame(table(all.ext.spp$region))
  if (nrow(all.reg.ext.spp)==0) { all.reg.ext.spp = data.frame(Var1=1:10, Freq=rep(NA,10))}
  
  reg.rich = merge(all.reg.rich, cl.reg.rich, by='Var1',all.x=T)
  names(reg.rich) = c('region','total.rich','clade.rich')
  
  reg.ext1 = merge(all.reg.ext.pops, cl.reg.ext.pops, by='Var1', all.x=T)
  names(reg.ext1) = c('Var1','total.ext.pops','clade.ext.pops')
  reg.ext2 = merge(reg.ext1, all.reg.ext.spp, by = 'Var1', all.x=T)
  reg.ext3 = merge(reg.ext2, cl.reg.ext.spp, by = 'Var1', all.x=T)
  names(reg.ext3)[c(1,4,5)] = c('region','total.ext.spp','clade.ext.spp')
  
  reg.rich.ext = merge(reg.rich, reg.ext3, by = 'region', all.x=T)
  reg.rich.ext$region = as.numeric(as.character(reg.rich.ext$region))
  reg.rich.ext = reg.rich.ext[order(reg.rich.ext$region),]
    
  reg.rich.thru.time = rbind(reg.rich.thru.time, cbind(time=rep(t,nrow(reg.rich.ext)), reg.rich.ext))
}

cols = rainbow(length(timeslices)+1)

pdf(paste(newplot_dir,"/sim",sim,"_focalclade",clade,"_thru_time_",Sys.Date(),".pdf",sep=""), height = 6, width = 10)

# Richness plots
par(mfrow=c(1,2),oma=c(0,0,2.5,0),mar=c(4,4,1,1))
plot(c(1,11),c(0,log10(max(reg.rich.thru.time$total.rich, na.rm=T))),type="n",xlab="Latitude",ylab="log10 Species richness")
sapply(timeslices,function(x) points(11 - reg.rich.thru.time$region[reg.rich.thru.time$time==x], 
                                     log10(reg.rich.thru.time$total.rich[reg.rich.thru.time$time==x]),
                                     type='l', lwd = 2, col=cols[which(timeslices==x)]))
par(new=T)
plot(c(1,11),c(0,log10(max(reg.rich.thru.time$total.rich, na.rm=T))),type="n",xlab="",ylab="",yaxt="n")
sapply(timeslices,function(x) points(11 - reg.rich.thru.time$region[reg.rich.thru.time$time==x], 
                                     log10(reg.rich.thru.time$clade.rich[reg.rich.thru.time$time==x]),
                                     type='l', lty='dotdash', lwd=4, col=cols[which(timeslices==x)]))
points(11 - clade.3011.pops$region[clade.3011.pops$time.of.origin==min(clade.3011.pops$time.of.origin)],
       0, pch=17, cex = 2, col = 'red')
legend('bottomleft',legend=timeslices,lty='solid',col=cols[1:length(timeslices)])
mtext("Species richness" , 3)

# Extinction plots
#par(mfrow=c(1,1),oma=c(0,0,2.5,0),mar=c(4,4,1,1))
plot(c(1,11),c(0,log10(max(reg.rich.thru.time$total.ext.pops, na.rm=T))),type="n",xlab="Latitude",ylab="log10 Extinct Populations")
sapply(timeslices,function(x) points(11 - reg.rich.thru.time$region[reg.rich.thru.time$time==x], 
                                     log10(reg.rich.thru.time$total.ext.pops[reg.rich.thru.time$time==x]),
                                     type='l', lwd = 2, col=cols[which(timeslices==x)]))
par(new=T)
plot(c(1,11),c(0,log10(max(reg.rich.thru.time$total.ext.pops, na.rm=T))),type="n",xlab="",ylab="",yaxt="n")
sapply(timeslices,function(x) points(11 - reg.rich.thru.time$region[reg.rich.thru.time$time==x], 
                                     log10(reg.rich.thru.time$clade.ext.pops[reg.rich.thru.time$time==x]),
                                     type='b', lty='dotdash', lwd=4, col=cols[which(timeslices==x)], cex=.3))

legend('topright',legend=c('All species','Focal clade'),lty=c('solid','dotdash'),lwd=c(2,4))
mtext("Population extinctions", 3)
mtext(paste("Sim = ",sim,"; focal cladeID = ",clade, sep=""), 3, outer=T, cex = 1.5)
dev.off()

