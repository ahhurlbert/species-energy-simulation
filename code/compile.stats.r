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

# below here the code is combining stats files generated using the analysis workflow on 10 Oct 2013, which includes a log distribution
# for the tree imbalance (i.e. beta) metric, hence the stats output file names
# energy gradient simulations were run to 100,000 time steps. Time scenario simulations effectively run to the 10,000 spp. richness limit
# included are 100 replicates for each scenario with each region of origin, so writing out 4 files

setwd("//constance/people/steg815/senc.analysis");

sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv"); head(sim.matrix);

Ttrop.sims = c(3465:3564); length(Ttrop.sims); # 100
Ttemp.sims = c(3565:3664); length(Ttemp.sims); # 100
Ktrop.sims = c(4065:4074,4185:4274); length(Ktrop.sims); # 100
Ktemp.sims = c(4075:4084,4275:4364); length(Ktemp.sims); # 100
all.sims = c(Ttrop.sims,Ttemp.sims,Ktrop.sims,Ktemp.sims); length(unique(all.sims)); # 400

stats.comp.Ttrop = numeric();
stats.comp.Ttemp = numeric();
stats.comp.Ktrop = numeric();
stats.comp.Ktemp = numeric();

for (i in all.sims) {
  
  if (sim.matrix$reg.of.origin[sim.matrix$sim.id == i] == 'tropical' & 
        sim.matrix$carry.cap[sim.matrix$sim.id == i] == 'off' &
        sim.matrix$energy.gradient[sim.matrix$sim.id == i] == 'off') {stats.comp.Ttrop = rbind(stats.comp.Ttrop,read.csv(paste("lbeta_Stats_sim",i,"_endtime_all_clades.csv",sep=""),header=T))};
  
  if (sim.matrix$reg.of.origin[sim.matrix$sim.id == i] == 'temperate' & 
        sim.matrix$carry.cap[sim.matrix$sim.id == i] == 'off' &
        sim.matrix$energy.gradient[sim.matrix$sim.id == i] == 'off') {stats.comp.Ttemp = rbind(stats.comp.Ttemp,read.csv(paste("lbeta_Stats_sim",i,"_endtime_all_clades.csv",sep=""),header=T))};
  
  if (sim.matrix$reg.of.origin[sim.matrix$sim.id == i] == 'tropical' & 
        sim.matrix$carry.cap[sim.matrix$sim.id == i] == 'on' &
        sim.matrix$energy.gradient[sim.matrix$sim.id == i] == 'on' &
        sim.matrix$max.time[sim.matrix$sim.id == i] == 1e+05) {stats.comp.Ktrop = rbind(stats.comp.Ktrop,read.csv(paste("lbeta_Stats_sim",i,"_endtime_all_clades.csv",sep=""),header=T))};
  
  if (sim.matrix$reg.of.origin[sim.matrix$sim.id == i] == 'temperate' & 
        sim.matrix$carry.cap[sim.matrix$sim.id == i] == 'on' &
        sim.matrix$energy.gradient[sim.matrix$sim.id == i] == 'on' &
        sim.matrix$max.time[sim.matrix$sim.id == i] == 1e+05) {stats.comp.Ktemp = rbind(stats.comp.Ktemp,read.csv(paste("lbeta_Stats_sim",i,"_endtime_all_clades.csv",sep=""),header=T))};
 
  print(i);
  
}; # end all.sims loop

head(stats.comp.Ttrop); dim(stats.comp.Ttrop); length(unique(stats.comp.Ttrop$sim)); # 100
head(stats.comp.Ttemp); dim(stats.comp.Ttemp); length(unique(stats.comp.Ttemp$sim)); # 100
head(stats.comp.Ktrop); dim(stats.comp.Ktrop); length(unique(stats.comp.Ktrop$sim)); # 100
head(stats.comp.Ktemp); dim(stats.comp.Ktemp); length(unique(stats.comp.Ktemp$sim)); # 100

write.csv(stats.comp.Ttrop,"SENC_Stats_100K_endtime_T.trop.csv",row.names=F,quote=F);
write.csv(stats.comp.Ttemp,"SENC_Stats_100K_endtime_T.temp.csv",row.names=F,quote=F);
write.csv(stats.comp.Ktrop,"SENC_Stats_100K_endtime_K.trop.csv",row.names=F,quote=F);
write.csv(stats.comp.Ktemp,"SENC_Stats_100K_endtime_K.temp.csv",row.names=F,quote=F);



