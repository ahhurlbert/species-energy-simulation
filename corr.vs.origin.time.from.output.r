analysis_dir = 'C:/SENCoutput/senc_analysis_130109'

sim.matrix = as.data.frame(read.csv("C:/Documents and Settings/Hurlbert/species-energy-simulation/SENC_Master_Simulation_Matrix.csv",header=T));
sim.matrix = subset(sim.matrix, !(energy.gradient=="on" & carry.cap=="off"))

which.sims = unique(sim.matrix$sim.id)

for (sim in which.sims) {
  
  output = read.csv(paste(analysis_dir,"/SENC_Stats_sim",sim,".csv",sep=""), header=T)
  ####################################################
  # Simulation summary plots
  ####################################################

  # clade.origin.corr.plot only if there are some output rows 
  if(!is.null(nrow(output))) {
    clade.origin.corr.plot(output, params.out, output.dir = analysis_dir)
  }
}