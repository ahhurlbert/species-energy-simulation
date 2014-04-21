setwd('//bioark.bio.unc.edu/hurlbertallen/manuscripts/frontierstropicaldiversity')
source('//bioark.bio.unc.edu/hurlbertallen/git/species-energy-simulation/code/supplemental_analysis_functions.r')

sim.matrix = read.csv('SENC_Master_Simulation_Matrix.csv')


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

# Function to plot regional abundance distributions (histograms) at the end of the simulation
# --sim: simulation ID
# --sim_dir: path in which sim output folders are stored
# --sim.matrix

sim.abun.dist = function(sim, sim_dir, sim.matrix) {
  simdata = output.unzip(sim_dir, sim)
  all.pops = simdata$all.populations
  params = sim.matrix[sim.matrix$sim.id == sim,]
  w = params$w
  carry.cap = params$carry.cap
  time = params$max.time
  
  pops.extant = subset(all.pops, time.of.extinction > time, 
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

# Probability that no species go extinct over a fixed period of time given that all species have 
# equal abundances under the overall constraint. (A more even abundance distribution results in more extinction.)
#p.no.extinct = function(gamma = 0.1, K, S, alpha = 1e-6) { ((1 - exp(-gamma*(K/S)))^S)^(1/(K*alpha)) }
#plot(S, p.no.extinct(K = K, S= S), type = 'l', ylab = 'Prob of no extinction')

# Probability that a species DOES go extinct prior to the appearance of the next new species.
# Expected waiting time until new species arises = 1/(K*alpha).
p.extinct = function(gamma = 0.1, K, S, alpha = 1e-6) { 1 - ((1 - exp(-gamma*(K/S)))^S)^(1/(K*alpha)) }

# This is a steep sigmoidal function as a function of S where the inflection point varies
# with gamma. Basically goes from 0 probability of any extinction for low S to perfect certainty
# of extinction. Alpha (speciation parameter) does not factor in, should only determine how 
# quickly this equilibrium value is achieved.

#At equilibrium under zero sum energy gradient scenario:

Ks = seq(4000, 40000, length.out = 12)[2:11] #because bins 0 and 11 were excluded
pinkred = colorRampPalette(c('pink','red','darkred','black'))
blues = colorRampPalette(c('lightblue', 'blue','darkblue'))

# In reds, the effect of varying K on equilibrial S is shown for alpha = 1e-6.
# In blues, it is shown for alpha = 1e-5
S = 1:600

pdf('effectsOf_AlphaGammaK_on_equilibRichness.pdf', height = 6, width = 18)
par(mfrow = c(1, 3), mar = c(6,6,1,1), cex.lab = 2, cex.axis = 2, las = 1, mgp = c(4, 1, 0))
plot(S, p.extinct(K = K, S= S), type = 'n', ylab = 'Prob of extinction before next new species',
     xlab = "")
sapply(1:10, function(x) points(S, p.extinct(K = Ks[x], S = S, alpha = 1e-6), 
                                type = 'l', col = pinkred(12)[x], lwd = 3))
sapply(1:10, function(x) points(S, p.extinct(K = Ks[x], S = S, alpha = 1e-5), 
                                type = 'l', col = blues(10)[x], lwd = 2, lty = 'dashed'))
abline(h= .5)
points(1:5*100, rep(.5, 5), pch = 15)
legend('bottomright', lwd = c(4,2,4,4), col = c('red','blue','gray90','gray40'), lty = c('solid','dashed','solid','solid'),
       legend = c('alpha = 1e-6', 'alpha = 1e-5', 'K = 7273', 'K = 36727'))

# The effect of varying alpha on equilibrial S is shown for K = 20000
alphas = seq(1e-6, 1e-4, length.out = 12)[2:11]
plot(S, p.extinct(K = 20000, S= S, alpha = 1e-5), type = 'n', xlab = 'Species richness', ylab = '')
sapply(1:10, function(x) points(S, p.extinct(K = 20000, S = S, alpha = alphas[x]), 
                                type = 'l', col = pinkred(12)[x], lwd = 3))
abline(h= .5)
points(1:5*100, rep(.5, 5), pch = 15)
legend('bottomright', lwd = 4, col = c('pink', 'darkred'), legend = c('alpha = 1e-5', 'alpha = 9e-5'))

# The effect of varying alpha on equilibrial S is shown for K = 20000
gammas = c(0.02, 0.05, 0.1, 0.2, 0.3)
plot(S, p.extinct(K = 20000, S= S, alpha = 1e-6), type = 'n', xlab = '', ylab = '')
sapply(1:5, function(x) points(S, p.extinct(K = 20000, S = S, alpha = 1e-6, gamma = gammas[x]), 
                                type = 'l', col = pinkred(6)[x], lwd = 3))
abline(h= .5)
points(1:5*100, rep(.5, 5), pch = 15)
legend('bottomright', lwd = 4, col = pinkred(6)[c(1,3,5)], 
       legend = c('gamma = 0.02', 'gamma = 0.1', 'gamma = 0.3'))
dev.off()
