setwd('//bioark.bio.unc.edu/hurlbertallen/manuscripts/frontierstropicaldiversity')
source('C:/species-energy-simulation/code/supplemental_analysis_functions.r')

sim.matrix = read.csv('C:/species-energy-simulation/SENC_Master_Simulation_Matrix.csv')


# function for adding region-specific alpha and carrying capacity to all.populations
alpha.K = function(sim.matrix, sim) {
  alpha = sim.matrix$alpha[sim.matrix$sim.id == sim]  #per capita speciation rate
  carry.cap = sim.matrix$carry.cap[sim.matrix$sim.id == sim]  #'on' if limits exist to community abundance; else 'off'
  energy.gradient = sim.matrix$energy.gradient[sim.matrix$sim.id == sim]  #'on' if maximum abundance varies across the gradient; else 'off'
  specn.gradient = sim.matrix$specn.gradient[sim.matrix$sim.id == sim] #'on' if specn rate is higher in tropics; 'off' if does not vary across gradient
  specn.factor = sim.matrix$specn.factor[sim.matrix$sim.id == sim] #factor by which tropical speciation rate is greater than temperate
  max.K = sim.matrix$max.K[sim.matrix$sim.id == sim]  #num of individuals that can be supported in region with the highest carrying capacity
  min.K = max.K/10  #num of individuals that can be supported in region with the lowest carrying capacity when there's a gradient
  num.of.bins = sim.matrix$num.of.bins[sim.matrix$sim.id == sim]  #number of spatial bins
  
  reg.alpha.K = data.frame(region = 0:num.of.bins)
  if (energy.gradient == 'on') { 
    reg.alpha.K$carry.cap = seq(min.K, max.K, length = num.of.bins + 1)
  } 
  if (energy.gradient == 'off') { 
    reg.alpha.K$carry.cap = rep(max.K, num.of.bins + 1)
  } 
  # regional speciation rates depending on whether they vary across the gradient
  if (specn.gradient == 'on') {
    reg.alpha.K$reg.alpha = seq(alpha/(specn.factor ^ .5), alpha*(specn.factor ^ .5), length = num.of.bins + 1)
  }
  if (specn.gradient == 'off') {
    reg.alpha.K$reg.alpha = rep(alpha, num.of.bins + 1)
  }
  return(reg.alpha.K)
}

sim.abun.dist = function(sim, sim_dir, sim.matrix) {
  simdata = output.unzip(sim_dir, sim)
  all.pops = simdata$all.populations
  params = sim.matrix[sim.matrix$sim.id == sim,]
  w = params$w
  carry.cap = params$carry.cap
  
  pops.extant = subset(all.pops, time.of.extinction > 100000, 
                       select = c('spp.name', 'region', 'env.opt','reg.env','carry.cap','time.of.sp.origin','time.of.origin'))
  pops.extant$pop.size = exp(-(pops.extant$env.opt - pops.extant$reg.env)^2/(w^2)) * 
    pops.extant$carry.cap
  reg.pop.sizes = as.data.frame(tapply(pops.extant$pop.size, pops.extant$region, FUN = sum))
  reg.pop.sizes = cbind(rownames(reg.pop.sizes), reg.pop.sizes) 
  colnames(reg.pop.sizes) = c("region", "reg.individuals")
  pops.extant = merge(pops.extant, reg.pop.sizes, by = "region")
  
  #The pop.multiplier assures that the sum of species abundances will not exceed carrying capacity
  pops.extant$pop.multiplier = pops.extant$carry.cap / pops.extant$reg.individuals
  pops.extant$pop.multiplier[pops.extant$pop.multiplier > 1] = 1
  
  # if there is no carrying capacity, then population size is not influenced by the actual 
  # community-level number of individuals, and so the sum across species is not modified (pop.multiplier = 1)
  if (carry.cap == 'off') {pops.extant$pop.multiplier = 1} 
  
  pops.extant$real.pop.size = round( pops.extant$pop.size * pops.extant$pop.multiplier , digits = 0)
  
  par(mfrow = c(3,4), mar = c(3,3,1,1), oma = c(0, 0, 2, 0))
  for (i in 1:10) { 
    hist(pops.extant$real.pop.size[pops.extant$region==i], xlab = "", main="", ylab="")
    mtext(paste("Sim", sim, "abundance distribution at end"), 3, outer = T)
  }
  return(pops.extant)
}

##################
time4565 = read.csv('raw_sim_output/sim4565_out/SENC_time.rich_sim4565.csv')
pops4565 = read.csv('raw_sim_output/sim4565_out/SENC_all.pops_sim4565.csv')
time5245 = read.csv('raw_sim_output/sim5245_out/SENC_time.rich_sim5245.csv')
pops5245 = read.csv('raw_sim_output/sim5245_out/SENC_all.pops_sim5245.csv')
time5445 = read.csv('raw_sim_output/sim5445_out/SENC_time.rich_sim5445.csv')
pops5445 = read.csv('raw_sim_output/sim5445_out/SENC_all.pops_sim5445.csv')
sim4065 = output.unzip('//bioark.bio.unc.edu/hurlbertallen/sencoutput/hurlbert_and_stegen_2014/raw_sim_output', 4065)
time4065 = sim4065$time.richness
pops4065 = sim4065$all.populations

data.out = c()
for (s in c(4065, 4565, 5245, 5445)) {
  temp = alpha.K(sim.matrix, s)
  finalrich = get(paste('time', s, sep = ''))
  temp$rich = finalrich$spp.rich[(nrow(finalrich) - 11):nrow(finalrich)]
  temp$sim = rep(s, 12)
  data.out = rbind(data.out, temp)
}

data.out = data.out[!data.out$region %in% c(0,11), ]
data.out$col = c(rep('red', 10), rep('lightblue', 10), rep('pink', 10), rep('lightgreen', 10))

par(mfrow = c(2, 1), mar = c(4,4,1,1))
plot(jitter(data.out$carry.cap), data.out$rich, pch = 16, xlab = "K",
     ylab = "S", col = data.out$col)
data.K22 = subset(data.out, carry.cap == 22000)
plot(data.K22$reg.alpha, data.K22$rich, pch = 16, xlab = "alpha", 
     ylab = "S", col = data.K22$col)

abun4065 = sim.abun.dist(sim.matrix = sim.matrix, sim_dir = '//bioark.bio.unc.edu/hurlbertallen/sencoutput/hurlbert_and_stegen_2014/raw_sim_output', sim = 4065)
par(mfrow=c(1,1), mar=c(4,4,1,1))
plot(log10(100001 - abun4065$time.of.origin), abun4065$real.pop.size, xlab = "Species age", ylab = "Abundance")
