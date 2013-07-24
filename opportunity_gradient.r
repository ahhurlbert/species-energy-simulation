source('unzipping_files.r')

sim_dir = "C:/SENCoutput"

#Sim with energy gradient and zero sum constraint and tropical origin
sim=3325
stats.3325 = read.csv(paste(sim_dir,'/SENC_stats_sim',sim,'.csv',sep=''), header=T)

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

#Identify clades that originated within 20% of 6000 timesteps and that came to span >=5 regions
f=.2
stats.3325.6k = subset(stats.3325, clade.origin.time > 6000*(1-f) & clade.origin.time < (1+f)*6000 & n.regions >=5)
#Clade 3193 has a strong reverse latitudinal richness gradient (r = 0.89) [r.env.rich = -0.89]
c = 3193
timeoforigin.6k = stats.3325.6k$clade.origin.time[stats.3325.6k$clade.id==c]
sub.clade = clade.members(c, sub.phylo, tip.labels=T)
subset.populations = subset(all.populations, spp.name %in% as.numeric(sub.clade));
#sub.populations is the subset of populations specific to a particular clade and timeslice
sub.populations = subset(subset.populations, time.of.origin <= t & time.of.extinction > t)
#Region of origin of clade
reg.of.origin = sub.populations[sub.populations$time.of.sp.origin == min(sub.populations$time.of.sp.origin), 'region']
clade.latrich = data.frame(table(sub.populations$region))
names(clade.latrich) = c('region','richness')
clade.latrich$region = as.numeric(as.character(clade.latrich$region))

# Function for characterizing the "opportunity gradient" at any point in time
oppgrad = function(t) {
  total.latgrad = subset(time.richness, time==30000, select = c('region','spp.rich'))
  timeoforigin.latgrad = subset(time.richness, time == t, select = c('region','spp.rich'))
  opportunity = merge(total.latgrad, timeoforigin.latgrad, by = 'region')
  names(opportunity)[2:3] = c('eq.rich','timeoforigin.rich')
  opportunity$opp = 1 - opportunity$timeoforigin.rich / opportunity$eq.rich
  opportunity = subset(opportunity, region > 0 & region < 11)
#  opp.grad = -cor(opportunity$region, opportunity$opp)
#  return(opp.grad)
}


#Identify clades with weak lat richness gradient and that came to span >=5 regions
w = .2
stats.3325.weak = subset(stats.3325, abs(r.env.rich) < w & n.regions >=5)
#Clade 2457
c = 2457
timeoforigin.weak = stats.3325.weak$clade.origin.time[stats.3325.weak$clade.id==c]

#Plot
opportunity = oppgrad(timeoforigin.weak)
par(mgp = c(3, 1, 0), mar = c(4, 4, 1, 4), oma = c(1, 1, 1, 1))
plot(11-opportunity$region,log10(opportunity$eq.rich), type='l', xlim=c(1,10), ylim = c(0,3),
     xlab = 'Latitudinal bin', ylab = 'log10 Species richness')
points(11-opportunity$region, log10(opportunity$timeoforigin.rich), type='l', col='red')
points(11-clade.latrich$region, log10(clade.latrich$richness), type = 'l', col='blue')
arrows(11-reg.of.origin, 0.25*max(log10(time.richness$spp.rich)), 11-reg.of.origin, 0)
par(new=T)
plot(11-opportunity$region, opportunity$opp, type = 'l', xlim = c(1,10), col='green', yaxt="n", xlab="", ylab="")
axis(4)
mtext("Proportion of equilibrial richness available", 4, line = 3)


# Relationship between "opportunity gradient" and latitudinal richness gradient for all subclades
# based on the "richness deficit" relative to equilibrial values at the time a subclade arose
times = unique(stats.3325$clade.origin.time)
opp.grads = sapply(times, function(x) {oppgrad(x))

opp.grad.times = data.frame(time = times, opp.grad = opp.grads)

stats.opp = merge(stats.3325, opp.grad.times, by.x='clade.origin.time', by.y='time',all.x=T)
plot(stats.opp$opp.grad, -stats.opp$r.env.rich) # No relationship between opportunity gradient and latitudinal gradient