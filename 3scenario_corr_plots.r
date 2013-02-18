# Manuscript figure
setwd('C:/Documents and Settings/Hurlbert/species-energy-simulation')


sim.matrix = read.csv('SENC_Master_Simulation_Matrix.csv',header=T)

source('time.env.corr.plots.r')
source('rep.plot.corrs.r')

# Pull out 6 subsets of sim.matrix: 3 treatments by 2 regions of origin
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

par(mfrow=c(3,4),mar=c(4,4,4,1),oma=c(3,1,1,1))
rep.plot(Kgrad.tr)
rep.plot(Kgrad.te)
rep.plot(dist.tr)
rep.plot(dist.te)
rep.plot(NoLim.tr)
rep.plot(NoLim.te)
