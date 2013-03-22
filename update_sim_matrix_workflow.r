# Workflow to update the master simulation matrix prior to running simulations

#(1) specify parameter space to explore
reg.of.origin = c('tropical','temperate'); # ancestor's region of origin, can be 'temperate' or 'tropical' or can vary between 0 and the 'num.of.bins', with 1 being temperate and the num.of.bins - 1 being tropical
w = c(3); # standard deviation of environmental niche function
alpha = c(0.000001); # per-individual mutation probability
alpha.to.beta = c(100); # muliplier (on alpha) to determine per-individual dispersal probability
sigma_E = c(1); # standard deviation of gaussian function from which mutant environmental optima are drawn
carry.cap = c('on'); # toggle to turn individual carrying capacity 'on' or 'off'
energy.gradient = c('off'); # toggle to turn a gradient in the number of individuals 'on' or 'off'
max.K = c(22000); # maximum individuals carry capacity in the region of highest abundance
num.of.bins = c(11); # one less than the actual number of bins, starting with bin #0. The two ends will be trimmed out at the end for analyses.
max.time = c(30000); # simulation time steps
max.richness = c(10000); # maximum total richness of the system
replicates = 1:100; # number of replicates for each param combination
disturb_frequency = c(100); # number of time steps between disturbance events
temperate_disturb_intensity = c(0.99); # fraction of individuals lost during a disturbance event in the temperate zone
tropical_disturb_intensity = c(0.75); # fraction of individuals lost during a disturbance event in the tropical zone

#(2) update master simulation matrix with chosen parameter combinations
setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation");
source('make_sim_matrix_fun.r');

sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T)); head(sim.matrix); dim(sim.matrix);
#sim.matrix$status = 'completed'; # comment in or out at needed
sim.matrix.new = make_sim_matrix(curr.matrix = sim.matrix,reg.of.origin=reg.of.origin,w=w,alpha=alpha,alpha.to.beta=alpha.to.beta,sigma_E=sigma_E,carry.cap=carry.cap,energy.gradient=energy.gradient,max.K=max.K,num.of.bins=num.of.bins,max.time=max.time,max.richness=max.richness,disturb_frequency=disturb_frequency,temperate_disturb_intensity=temperate_disturb_intensity,tropical_disturb_intensity=tropical_disturb_intensity,replicates=replicates)
head(sim.matrix.new[sim.matrix.new$status=='to.run',]); dim(sim.matrix.new[sim.matrix.new$status=='to.run',]);
write.csv(sim.matrix.new,"SENC_Master_Simulation_Matrix.csv",row.names=F,quote=F);
