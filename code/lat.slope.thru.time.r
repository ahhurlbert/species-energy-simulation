setwd('z:/git/species-energy-simulation/')
source('code/supplemental_analysis_functions.r')
sim3465 = output.unzip('z:/sencoutput/hurlbert_and_stegen_2014/raw_sim_output', 3465)
time3465 = sim3465$time.richness
time3465 = time3465[!time3465$region %in% c(0,11),]

sim4065 = output.unzip('z:/sencoutput/hurlbert_and_stegen_2014/raw_sim_output', 4065)
time4065 = sim4065$time.richness
time4065 = time4065[!time4065$region %in% c(0,11),]

sim5525 = output.unzip('z:/manuscripts/frontierstropicaldiversity/raw_sim_output', 5525)
time5525 = sim5525$time.richness
time5525 = time5525[!time5525$region %in% c(0,11),]

sim5625 = output.unzip('z:/manuscripts/frontierstropicaldiversity/raw_sim_output', 5625)
time5625 = sim5625$time.richness
time5625 = time5625[!time5625$region %in% c(0,11),]


# function for calculating the ratio of richness in the tropics to richness
# in the temperate most region that is occupied
trop.temp.ratio = function(temp.time.rich) {
  ratio = temp.time.rich$spp.rich[temp.time.rich$region == 10] / 
    temp.time.rich$spp.rich[temp.time.rich$region == min(temp.time.rich$region[temp.time.rich$spp.rich != 0])]
}

times = 1:1000*100
output = c()
for (i in times) {
  temp4065 = subset(time4065, time == i)
  temp5525 = subset(time5525, time == i)
  temp5625 = subset(time5625, time == i)
  lm4065 = lm(spp.rich~region, data = temp4065)
  lm5525 = lm(spp.rich~region, data = temp5525)
  lm5625 = lm(spp.rich~region, data = temp5625)
  ratio4065 = trop.temp.ratio(temp4065)
  ratio5525 = trop.temp.ratio(temp5525)
  ratio5625 = trop.temp.ratio(temp5625)
  output = rbind(output, c(-lm4065$coefficients[2],
                           -lm5525$coefficients[2],
                           -lm5625$coefficients[2],
                           ratio4065,
                           ratio5525,
                           ratio5625))
}
latslope = data.frame(cbind(times, output))
names(latslope) = c('time', 'slope4065', 'slope5525', 'slope5625', 'ratio4065', 'ratio5525', 'ratio5625')


times2 = 1:100*2
out = c()
for (i in times2) {
  temp3465 = subset(time3465, time == i)
  lm3465 = lm(spp.rich~region, data = temp3465)
  out = c(out, -lm3465$coefficients[2])
}
slope3465 = data.frame(cbind(times2, out))
names(slope3465) = c('time', 'slope')

pdf('z:/manuscripts/frontierstropicaldiversity/lat_slope_thru_time_4scenarios.pdf', height = 5, width = 6)
par(mar = c(4,5,1,1), las = 1, cex.lab = 1.5)
plot(slope4065 ~ time, data = latslope, type = 'l', xlab = 'Time', ylab = 'Latitudinal gradient slope', ylim = c(-53,0))
points(latslope$time, latslope$slope5525, type = 'l', col = 'red')
points(latslope$time, latslope$slope5625, type = 'l', col = 'blue')
points(300*slope3465$time, slope3465$slope, type = 'l', col = 'springgreen1')
text.cex = 1
text(6e4, -3, "Speciation gradient", col = 'red', cex = text.cex)
text(6e4, -19, "Disturbance gradient", col = 'blue', cex = text.cex)
text(6e4, -42, "Energy gradient", cex = text.cex)
text(6e4, -51, "Pure niche\nconservatism", col = 'springgreen1', cex = text.cex)
dev.off()

if(0){
plot(log10(ratio4065) ~ time, data = latslope, type = 'l', xlab = 'Time', ylab = 'Tropical / temperate richness', ylim = c(0, 2.5))
points(latslope$time, log10(latslope$ratio5525), type = 'l', col = 'red')
points(latslope$time, log10(latslope$ratio5625), type = 'l', col = 'blue')
text(6e4, -3, "Speciation gradient", col = 'red')
text(6e4, -19, "Disturbance gradient", col = 'blue')
text(6e4, -42, "Energy gradient")
}