#!/usr/bin/env Rscript

#sim = commandArgs();
#sim = as.numeric(sim[length(sim)]);

Allen = 0;
partial.analysis = 1; # toggle to determine whether we're looking at all sims or just some

#New parameter for taking into account which of us is running this code
if(Allen==1) {
  setwd('c:/documents and settings/hurlbert/species-energy-simulation')
  Rlib.location = "C:/program files/R/R-2.15.2/library"
  sim_dir = "C:/SENCoutput"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses"
} else {
  setwd('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation')
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever all of your zipped output files are
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever you want to store the results of these analyses
}

# Simulation workflow

#(2) load analysis functions

source('pops.through.time.r');
source('unzipping_files.r');


#(3) read in master simulation matrix with chosen parameter combinations;
sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));

#(4) start analyses based on value of 'sim' which draws parameter values from sim.matrix
if (partial.analysis == 0) {which.sims = 1:max(sim.matrix$sim.id)};
if (partial.analysis == 1) {which.sims = c(3465,3665,3865)}; # which.sims = c(read.csv(paste(analysis_dir,"/sims.to.analyze.csv",sep=""))$x)

for (sim in which.sims) {
  
  # (5) read in simulation results for specified simulation from the output zip file
  sim.results = output.unzip(sim_dir,sim)
  
  all.populations = sim.results$all.populations
  head(all.populations)
  
  max.time = max(all.populations$time.of.extinction[all.populations$time.of.extinction < 30001]);
  
  num.of.samples = 200;
  
  pops.comp = array(0,dim=c(10,5,num.of.samples));
  pops.comp[,1,] = 1:10;
  
  counter = 1;
  
  for (t in seq(1,max.time,length=num.of.samples)) {
    
    pops.temp = pops.through.time(all.populations,t); head(pops.temp);
    pops.temp = merge(pops.comp[,1,counter],pops.temp,by=1,all=T)
    pops.comp[1:nrow(pops.temp),,counter] = as.matrix(pops.temp);
    rm('pops.temp');
    
    counter = counter + 1;
    print(c(t,counter))
  }
  
  # names of columns in array: region extinct.pops total.pops time.in.region.pops extant.pops
  
  sim.matrix[sim.matrix$sim.id==sim,];
  
  pdf(paste(analysis_dir,'/Pops.Thru.Time_',sim,'.pdf',sep=""),height=10,width=5);
    par(mfrow=c(3,1),pty="s");
  
    bin = 10; col.id = 3; plot(log10(pops.comp[bin,col.id,]) ~ pops.comp[bin,4,],ylab="Log10 (Cumulative Populations)",xlab="Time",col=2,typ="l",lwd=2,xlim=c(0,max.time),ylim=c(0,max(log10(pops.comp[,col.id,]),na.rm=T)));
    bin = 1; if(sim == 3465) {bin = 7}; col.id = 3; points(log10(pops.comp[bin,col.id,]) ~ pops.comp[bin,4,],col=4,typ="l",lwd=2);
  
    bin = 10; col.id = 2; plot(log10(pops.comp[bin,col.id,]) ~ pops.comp[bin,4,],ylab="Log10 (Cumulative Extinctions)",xlab="Time",col=2,typ="l",lwd=2,xlim=c(0,max.time),ylim=c(0,max(log10(pops.comp[,col.id,]),na.rm=T)));
    bin = 1; if(sim == 3465) {bin = 7}; points(log10(pops.comp[bin,col.id,]) ~ pops.comp[bin,4,],col=4,typ="l",lwd=2);
  
    bin = 10; col.id = 5; plot(log10(pops.comp[bin,col.id,]) ~ pops.comp[bin,4,],ylab="Log10 (Species Richness)",xlab="Time",col=2,typ="l",lwd=2,xlim=c(0,max.time),ylim=c(0,max(log10(pops.comp[,col.id,]),na.rm=T)));
    bin = 1; if(sim == 3465) {bin = 7}; points(log10(pops.comp[bin,col.id,]) ~ pops.comp[bin,4,],col=4,typ="l",lwd=2);
  
  dev.off()
  
}

