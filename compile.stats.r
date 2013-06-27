# just rbinding stats output into a file for a given scenario and a given region of origin
# currently written for time slice using the energy gradient sims only

setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204");

sim.matrix = read.csv("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation/SENC_Master_Simulation_Matrix.csv"); head(sim.matrix);

time.slice = 5459;

foo = list.files(pattern='SENC_Stats_sim'); head(foo); length(foo);
foo2 = foo[grep('time5459',foo)]; head(foo2);
foo3 = sub('SENC_Stats_sim','',foo2); head(foo3);
K.sims = as.numeric(sub('_time5459.csv','',foo3)); head(K.sims); length(K.sims);


stats.comp.trop = numeric();
stats.comp.temp = numeric();

for (i in K.sims) {
  
  if (sim.matrix$reg.of.origin[sim.matrix$sim.id == i] == 'tropical') {stats.comp.trop = rbind(stats.comp.trop,read.csv(paste("SENC_Stats_sim",i,"_time",time.slice,".csv",sep=""),header=T))};
  if (sim.matrix$reg.of.origin[sim.matrix$sim.id == i] == 'temperate') {stats.comp.temp = rbind(stats.comp.temp,read.csv(paste("SENC_Stats_sim",i,"_time",time.slice,".csv",sep=""),header=T))};
  
}

head(stats.comp.trop); dim(stats.comp.trop); length(unique(stats.comp.trop$sim));
head(stats.comp.temp); dim(stats.comp.temp); length(unique(stats.comp.temp$sim));

write.csv(stats.comp.trop,"SENC_Stats_K.slice.sims.trop.csv",row.names=F,quote=F);
write.csv(stats.comp.temp,"SENC_Stats_K.slice.sims.temp.csv",row.names=F,quote=F);