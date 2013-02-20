# Plot boxplots of gamma (phylogenetic tree shape statistic) across 3 scenarios
# (K gradient, disturbance gradient, or no limit) and 2 regions of origin.
sim_dir = "C:/SENCoutput"
analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/"
repo_dir = "C:/Documents and Settings/Hurlbert/species-energy-simulation"

sim.matrix = read.csv(paste(repo_dir,'/SENC_Master_Simulation_Matrix.csv',sep=""),header=T)

# Scenarios
dist.tr = subset(sim.matrix, max.K==22000 & reg.of.origin=='tropical')
dist.te = subset(sim.matrix, max.K==22000 & reg.of.origin=='temperate')
Kgrad.tr = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='tropical' & carry.cap=='on' & energy.gradient=='on' & max.richness==10000)
Kgrad.te = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='temperate' & carry.cap=='on' & energy.gradient=='on' & max.richness==10000)
NoLim.tr = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='tropical' & carry.cap=='off' & energy.gradient=='off' & max.richness==10000)
NoLim.te = subset(sim.matrix, w==3 & sigma_E==1 & max.K==40000 & alpha==1e-6 & beta==1e-4 & disturb_frequency==0
                  & reg.of.origin=='temperate' & carry.cap=='off' & energy.gradient=='off' & max.richness==10000)

# Function rbinds together latitudinal richness gradient data across specified simulations from the most recent timestep
recent.lat.grad = function(sim.ids,cut.ends=T) {
  latgrad.sim.reps = c()
  for (i in sim.ids) {
    tmp = read.csv(paste(sim_dir,'/SENC_time.rich_sim',i,'.csv',sep=''), header=T)
    max.time = max(tmp$time)
    tmp2 = tmp[tmp$time==max.time,]
    tmp2$max.time = max.time
    tmp2$sim = i
    latgrad.sim.reps = rbind(latgrad.sim.reps, tmp2)
  }  
  if(cut.ends) { latgrad.sim.reps = subset(latgrad.sim.reps, region %in% 1:10)} #eliminate regions 0 and 11 from output
  return(latgrad.sim.reps)
}

old.lat.grad = function(sim.ids, time, cut.ends = T) {
  latgrad.sim.reps = c()
  for (i in sim.ids) {
    tmp = read.csv(paste(sim_dir,'/SENC_time.rich_sim',i,'.csv',sep=''), header=T)
    tmp2 = tmp[tmp$time==time,]
    tmp2$time = time
    tmp2$sim = i
    latgrad.sim.reps = rbind(latgrad.sim.reps, tmp2)
  }  
  if(cut.ends) { latgrad.sim.reps = subset(latgrad.sim.reps, region %in% 1:10)} #eliminate regions 0 and 11 from output
  return(latgrad.sim.reps)
}

latgrad.dist.tr = recent.lat.grad(dist.tr$sim.id)
latgrad.dist.te = recent.lat.grad(dist.te$sim.id)
latgrad.Kgrad.tr = recent.lat.grad(Kgrad.tr$sim.id)
latgrad.Kgrad.te = recent.lat.grad(Kgrad.te$sim.id)
latgrad.NoLim.tr = recent.lat.grad(NoLim.tr$sim.id)
latgrad.NoLim.te = recent.lat.grad(NoLim.te$sim.id)

oldgrad.dist.tr = old.lat.grad(dist.tr$sim.id, time = 2000)
oldgrad.dist.te = old.lat.grad(dist.te$sim.id, time = 2000)
oldgrad.Kgrad.tr = old.lat.grad(Kgrad.tr$sim.id, time = 3000)
oldgrad.Kgrad.te = old.lat.grad(Kgrad.te$sim.id, time = 3000)
oldgrad.NoLim.tr = old.lat.grad(NoLim.tr$sim.id, time = 140)
oldgrad.NoLim.te = old.lat.grad(NoLim.te$sim.id, time = 140)

pdf(paste(analysis_dir,'/summaryplots/3scenarios_lat_grads.pdf',sep=''), height = 4, width = 8)
par(mfrow=c(1,3), mar = c(4,4,2,1), oma = c(2,1,1,1))
## Time scenario
plot(c(1,10),c(0,max(latgrad.NoLim.tr$spp.rich)),type="n",xlab="Latitude",ylab="Species richness", main="Time")
#recent gradient
sapply(NoLim.tr$sim.id, function(x)
       points(11 - latgrad.NoLim.tr[latgrad.NoLim.tr$sim==x,'region'], 
              latgrad.NoLim.tr[latgrad.NoLim.tr$sim==x,'spp.rich'], 
              type = 'l', lwd=2, col = 'gray80'))
sapply(NoLim.te$sim.id, function(x)
  points(11 - latgrad.NoLim.te[latgrad.NoLim.te$sim==x,'region'], 
         latgrad.NoLim.te[latgrad.NoLim.te$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'pink'))
#old gradient
sapply(NoLim.tr$sim.id, function(x)
  points(11 - oldgrad.NoLim.tr[oldgrad.NoLim.tr$sim==x,'region'], 
         oldgrad.NoLim.tr[oldgrad.NoLim.tr$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'gray40', lty='dashed'))
sapply(NoLim.te$sim.id, function(x)
  points(11 - oldgrad.NoLim.te[oldgrad.NoLim.te$sim==x,'region'], 
         oldgrad.NoLim.te[oldgrad.NoLim.te$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'indianred4', lty='dashed'))

## Disturbance scenario
plot(c(1,10),c(0,max(latgrad.dist.tr$spp.rich)),type="n",xlab="Latitude",
     ylab="Species richness", main="Disturbance")
#recent gradient
sapply(dist.tr$sim.id, function(x)
  points(11 - latgrad.dist.tr[latgrad.dist.tr$sim==x,'region'], 
         latgrad.dist.tr[latgrad.dist.tr$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'gray80'))
sapply(dist.te$sim.id, function(x)
  points(11 - latgrad.dist.te[latgrad.dist.te$sim==x,'region'], 
         latgrad.dist.te[latgrad.dist.te$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'pink'))

#old gradient
sapply(dist.tr$sim.id, function(x)
  points(11 - oldgrad.dist.tr[oldgrad.dist.tr$sim==x,'region'], 
         oldgrad.dist.tr[oldgrad.dist.tr$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'gray40', lty='dashed'))
sapply(dist.te$sim.id, function(x)
  points(11 - oldgrad.dist.te[oldgrad.dist.te$sim==x,'region'], 
         oldgrad.dist.te[oldgrad.dist.te$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'indianred4', lty='dashed'))

## K gradient scenario
plot(c(1,10),c(0,max(latgrad.Kgrad.tr$spp.rich)),type="n",xlab="Latitude",
     ylab="Species richness", main="Energetic constraints")
#recent gradient
sapply(Kgrad.tr$sim.id, function(x)
  points(11 - latgrad.Kgrad.tr[latgrad.Kgrad.tr$sim==x,'region'], 
         latgrad.Kgrad.tr[latgrad.Kgrad.tr$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'gray80'))
sapply(Kgrad.te$sim.id, function(x)
  points(11 - latgrad.Kgrad.te[latgrad.Kgrad.te$sim==x,'region'], 
         latgrad.Kgrad.te[latgrad.Kgrad.te$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'pink'))

#old gradient
sapply(Kgrad.tr$sim.id, function(x)
  points(11 - oldgrad.Kgrad.tr[oldgrad.Kgrad.tr$sim==x,'region'], 
         oldgrad.Kgrad.tr[oldgrad.Kgrad.tr$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'gray40', lty='dashed'))
sapply(Kgrad.te$sim.id, function(x)
  points(11 - oldgrad.Kgrad.te[oldgrad.Kgrad.te$sim==x,'region'], 
         oldgrad.Kgrad.te[oldgrad.Kgrad.te$sim==x,'spp.rich'], 
         type = 'l', lwd=2, col = 'indianred4', lty='dashed'))
dev.off()
