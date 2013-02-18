sim_dir = 'C:/SENCoutput'

rep.plot = function(sub.sim.matrix) {
  sim.uniq1 = unique(sub.sim.matrix[,c(3:9,15,16)])
  for (u in 1:nrow(sim.uniq1)) {
    sim.reps1 = subset(sub1, w == sim.uniq1$w[u] & sigma_E == sim.uniq1$sigma_E[u], select = sim.id)
    stats.output1 = c()
    for (r in 1:nrow(sim.reps1)) {
      temp = read.csv(paste(sim_dir,'/SENC_Stats_sim',sim.reps1[r,1],'.csv',sep=''),header=T)
      stats.output1 = rbind(stats.output1,temp)
    }
    clade.origin.corr.plot.simple(stats.output1, sim.matrix[sim.matrix$sim.id==sim.reps1[1,],], min.num.data.pts = 10, 
                                min.num.spp.per.clade = 30, min.num.regions = 5, output.dir )
  }
}

par(mfrow=c(2,3), mar = c(4,4,1,1), oma = c(3,1,3,1))
sub1 = subset(sim.matrix, w==3 & sigma_E ==1 & alpha==1e-6 & beta == 1e-4 & replicate>0 & disturb_frequency==0 & 
                carry.cap=="on" & energy.gradient=="on" & reg.of.origin=="tropical")
rep.plot(sub1)

sub2 = subset(sim.matrix, w==3 & sigma_E ==1 & alpha==1e-6 & beta == 1e-4 & replicate>0 & disturb_frequency==0 & 
                carry.cap=="on" & energy.gradient=="on" & reg.of.origin=="temperate")
rep.plot(sub2)
