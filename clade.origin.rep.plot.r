sim_dir = 'C:/SENCoutput/senc.out.130115'
output.dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/senc.out.130115'

sim.matrix = read.csv('C:/SENCoutput/senc.out.130115/sim.matrix.output_2013-01-16.csv',header=T)

#stats.output = compile.firstlines('C:/SENCoutput/senc.out.130115','SENC_Stats_sim')

w.sigma.layout.plot = function(reg, K, grad) {
  sub1 = subset(sim.matrix, reg.of.origin == reg & carry.cap == K & energy.gradient == grad)
  sim.uniq1 = unique(sub1[,3:9])
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



pdf(paste(output.dir,'/cladevscorr_10reps_',Sys.Date(),'.pdf',sep=''), height = 10, width = 8)
par(mfrow=c(4,3), oma=c(0,0,4,0), mar = c(3,4,1.5,1))
w.sigma.layout.plot('tropical','on','on')
w.sigma.layout.plot('tropical','off','off')
w.sigma.layout.plot('temperate','on','on')
w.sigma.layout.plot('temperate','off','off')
dev.off()
