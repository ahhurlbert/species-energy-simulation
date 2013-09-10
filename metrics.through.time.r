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
  sim_dir = "C:/SENCoutput"
  analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots"
} else {
  setwd('C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/species-energy-simulation')
  sim_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever all of your zipped output files are
  analysis_dir = "C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204" #wherever you want to store the results of these analyses
}

trop.sims = 4065:4067
temp.sims = 4075:4077

sim.matrix = read.csv("SENC_Master_Simulation_Matrix.csv",header=T);

trop.metrics = array(-999,c(length(trop.sims),100,28));
temp.metrics = array(-999,c(length(temp.sims),100,28));

#Need same number of sims in trop.sims and temp.sims
for (i in 1:length(trop.sims)) {
  
  sim1 = trop.sims[i];
  sim2 = temp.sims[i]
  temp1 = read.csv(paste(sim_dir,"/SENC_Stats_sim",sim1,"_mult_times.csv",sep=""),header=T)
  temp1$r.lat.rich = -temp1$r.env.rich; 
  trop.metrics[i,,] = as.matrix(temp1); 
  trop.metrics = array(trop.metrics, dim = dim(trop.metrics), dimnames = list(NULL, NULL, colnames(temp1)));
  temp2 = read.csv(paste(sim_dir,"/SENC_Stats_sim",sim2,"_mult_times.csv",sep=""),header=T);
  temp2$r.lat.rich = -temp2$r.env.rich;
  temp.metrics[i,,] = as.matrix(temp2);
  temp.metrics = array(temp.metrics, dim = dim(temp.metrics), dimnames = list(NULL, NULL, colnames(temp2)));
  rm(list = c('temp1', 'temp2'));
  
}

temp.metrics.mean = data.frame(apply(temp.metrics, 2:3, mean))
temp.metrics.sd = data.frame(apply(temp.metrics, 2:3, function(x) var(x)^.5))

trop.metrics.mean = data.frame(apply(trop.metrics, 2:3, mean))
trop.metrics.sd = data.frame(apply(trop.metrics, 2:3, function(x) var(x)^.5))

pdf(paste(analysis_dir,'/metrics_thru_time_',Sys.Date(),'.pdf',sep=""), height = 6, width = 8)
par(mfrow = c(2, 2), mar = c(3, 6, 1, 1), oma = c(5, 0, 0, 0), cex.lab = 2, las = 1, cex.axis = 1.3, mgp = c(4,1,0))
y.labels = c('Global richness', expression(italic(r)[latitude-richness]), expression(gamma), expression(italic(r)[MRD-richness]))
metric.names = c('global.richness','r.lat.rich', 'gamma.stat','r.MRD.rich')
for (j in 1:4) {
  curr.metric = metric.names[j]
  plot(trop.metrics.mean$time/1000, trop.metrics.mean[, curr.metric], xlim = c(0, max(trop.metrics.mean$time/1000)), 
       ylim = range(c(trop.metrics[, , curr.metric], temp.metrics[, , curr.metric])), type = "n",
       ylab = y.labels[j], xlab = "")
  polygon(c(trop.metrics.mean$time/1000, rev(trop.metrics.mean$time/1000)), 
          c(trop.metrics.mean[, curr.metric] - trop.metrics.sd[, curr.metric], 
            rev(trop.metrics.mean[, curr.metric] + trop.metrics.sd[, curr.metric])), 
          col = rgb(.8, 0, 0, .3), border = NA)
  polygon(c(temp.metrics.mean$time/1000, rev(temp.metrics.mean$time/1000)), 
          c(temp.metrics.mean[, curr.metric] - temp.metrics.sd[, curr.metric], 
            rev(temp.metrics.mean[, curr.metric] + temp.metrics.sd[, curr.metric])), 
          col = rgb(0, 0, .8, .3), border = NA)
  points(trop.metrics.mean$time/1000, trop.metrics.mean[, curr.metric], type = 'l', col = 'red', lwd = 3)
  points(temp.metrics.mean$time/1000, temp.metrics.mean[, curr.metric], type = 'l', col = 'blue', lwd = 3)
  if(curr.metric == 'gamma.stat') { abline(h = 0, lty = 'dashed')}
} 
mtext("Time (x1000)", 1, outer=T, cex = 1.75, line = 2)                                      
dev.off()
