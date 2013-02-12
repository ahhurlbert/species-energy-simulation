all.stats = compile.firstlines(sim_dir,"SENC_Stats_sim")
mama = merge(sim.matrix,all.stats,by.x='sim.id',by.y='sim',all.x=T)
mama2 = subset(mama, sim.id>2924)

pdf('//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses/summaryplots/disturbance_effect_on_gamma_Bk.pdf',
    height = 8, width = 6)
par(mfrow=c(2,1), mar = c(4,4,.5,1), oma=c(2,1,2,1))
boxplot(mama2$BK.env.y ~ mama2$disturb_frequency + mama2$temperate_disturb_intensity, xaxt="n", ylab = "Blomberg's K")
axis(1, labels = c(0, 0.75, 0.85, 0.95), at = c(1,7,11,15))
mtext(c("dist",0, 10,50,100, 10, 50, 100, 10, 50, 100), 1, at = c(-1, 1, 6:8, 10:12, 14:16), line = 2)
mtext("Sims 2925:3324, w = 3, sigma = 1, \ndisp = 1e-4, specn = 1e-6", outer=T)
boxplot(mama2$gamma.stat ~ mama2$disturb_frequency + mama2$temperate_disturb_intensity, xaxt = "n", ylab = "Gamma")
axis(1, labels = c(0, 0.75, 0.85, 0.95), at = c(1,7,11,15))
mtext(c("freq",0, 10,50,100, 10, 50, 100, 10, 50, 100), 1, at = c(-1, 1, 6:8, 10:12, 14:16), line = 2)
mtext("Disturbance regime", 1, outer=T)
abline(h=0,col='red',lty='dashed')
dev.off()

