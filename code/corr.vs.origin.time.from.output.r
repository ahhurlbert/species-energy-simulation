output_dir = 'C:/SENCoutput/senc_analysis_130109'
params_dir = 'C:/SENCoutput'
stats_dir = 'C:/SENCoutput/senc_analysis_121220'

sim.matrix = as.data.frame(read.csv("C:/Documents and Settings/Hurlbert/species-energy-simulation/SENC_Master_Simulation_Matrix.csv",header=T));
#subset out nonsensical parameter combos of energy gradient with no carrying capacity
sim.matrix = subset(sim.matrix, !(energy.gradient=="on" & carry.cap=="off"))

which.sims = unique(sim.matrix$sim.id[sim.matrix$sim.id!=0])

for (sim in which.sims) {
  
  output = read.csv(paste(stats_dir,"/SENC_Stats_sim",sim,".csv",sep=""), header=T)
  params.out = read.csv(paste(params_dir,"/SENC_params.out_sim",sim,".csv",sep=""),header=T)
  ####################################################
  # Simulation summary plots
  ####################################################

  # clade.origin.corr.plot only if there are some output rows 
  if(!is.null(nrow(output))) {
    clade.origin.corr.plot(output, params.out, output.dir = output_dir)
  }
}