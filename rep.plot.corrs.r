sim_dir = 'C:/SENCoutput'

rep.plot = function(sub.sim.matrix) {
  stats.output1 = c()
  for (r in 1:nrow(sub.sim.matrix)) {
    temp = read.csv(paste(sim_dir,'/SENC_Stats_sim',sub.sim.matrix$sim.id[r],'.csv',sep=''),header=T)
    stats.output1 = rbind(stats.output1,temp)
  }
  time.env.corr.plot(stats.output1, sim.params = sub.sim.matrix[1,], min.num.data.pts = 10, 
                                min.num.spp.per.clade = 30, min.num.regions = 5)
}

par(mfrow=c(2,3), mar = c(4,4,1,1), oma = c(3,1,3,1))
sub1 = subset(sim.matrix, w==3 & sigma_E ==1 & alpha==1e-6 & beta == 1e-4 & replicate>0 & disturb_frequency==0 & 
                carry.cap=="on" & energy.gradient=="on" & reg.of.origin=="tropical")
rep.plot(sub1)

sub2 = subset(sim.matrix, w==3 & sigma_E ==1 & alpha==1e-6 & beta == 1e-4 & replicate>0 & disturb_frequency==0 & 
                temperate_disturb_intensity==0 & carry.cap=="on" & energy.gradient=="on" & reg.of.origin=="temperate")
rep.plot(sub2)
