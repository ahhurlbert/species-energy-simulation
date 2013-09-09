## code for relating various metrics of the root clade to time and making plots

drop.extinct.tips = function(phy){
  tol<-1e-6
  temp<-diag(vcv(phy))
  pruned.phy<-drop.tip(phy,names(temp[temp<(max(temp)-tol)]))
  return(pruned.phy)
}

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


which.sims = c(4065,4075);

sim.matrix = as.data.frame(read.csv("SENC_Master_Simulation_Matrix.csv",header=T));

metrics = array(-999,c(length(which.sims),100,27));

for (i in 1:length(which.sims)) {
  
  sim = which.sims[i];
  temp.metrics = read.csv(paste(sim_dir,"/SENC_Stats_sim",sim,"_mult_times.csv",sep=""),header=T);
  metrics[i,,] = as.matrix(temp.metrics);
  metrics = array(metrics,dim=dim(metrics),dimnames = list(NULL,NULL,colnames(temp.metrics)));
  rm('temp.metrics');
  
}

pdf(paste(analysis_dir,'/senc.metrics.thru.time.pdf',sep=""))

for (curr.metric in c('global.richness','gamma.stat','r.time.rich','r.env.rich','r.MRD.rich','r.PSV.rich','r.env.MRD','r.env.PSV','r.ext.reg','r.rich.ext','BK.env','BK.reg')) {

  plot(metrics[1,,curr.metric] ~ metrics[1,,'time'],xlim=c(0,max(metrics[,,'time'])),ylim=range(metrics[,,curr.metric]),ylab=curr.metric,xlab='Time',typ='n');
  abline(v=30000,lty=2)

  for (i in 1:length(which.sims)) {
  
    sim = which.sims[i];
    if (sim.matrix$reg.of.origin[sim.matrix$sim.id == sim] == 'tropical') {col.use = 2};
    if (sim.matrix$reg.of.origin[sim.matrix$sim.id == sim] == 'temperate') {col.use = 4};
    points(metrics[i,,curr.metric] ~ metrics[1,,'time'],col=col.use);

  }
  
} # end curr.metric loop

dev.off()
