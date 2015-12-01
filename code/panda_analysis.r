# Code for analyzing the simulated phylogenies from 
# Hurlbert & Stegen 2014a,b that were generated under
# different eco-evolutionary dynamics using Morlon et al.'s
# (2010, PLoS Biology) RPANDA package

# Set working directory to species-energy-simulation repo

library(RPANDA)
library(geiger)

multi.panda.fit = function(tree) {
  
  # Model 1
  m1 = tryCatch(
    {
      fit_coal_cst(tree, tau0 = 1e-4, gamma = -1, cst.rate = TRUE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 1:", cond))
      list(model = "Equilibrium constant rate", LH = NA, aicc = NA, tau0 = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 1:", cond))
    },
    finally = message("Model 1 completed")
  )
  
  # Model 2  
  m2 = tryCatch(
    {
      fit_coal_cst(tree, tau0 = 1e-4, gamma = -1, cst.rate = FALSE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 2:", cond))
      list(model = "Equilibrium variable rate", LH = NA, aicc = NA, tau0 = NA, gamma = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 2:", cond))
    },
    finally = message("Model 2 completed")
  )
  
  # Model 3  
  m3 = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = TRUE, cst.mu = TRUE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 3:", cond))
      list(model = "Birth-death constant rates (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA, mu0 = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 3:", cond))
    },
    finally = message("Model 3 completed")
  )
  
  # Model 4a
  m4a = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = FALSE, cst.mu = TRUE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 4a:", cond))
      list(model = "Birth-death varying speciation constant extinction (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA, mu0 = NA, alpha = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 4a:", cond))
    },
    finally = message("Model 4a completed")
  )
  
  # Model 4b
  m4b = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = TRUE, cst.mu = FALSE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 4b:", cond))
      list(model = "Birth-death constant speciation varying extinction (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA, mu0 = NA, beta = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 4b:", cond))
    },
    finally = message("Model 4b completed")
  )

  # Model 4c
  m4c = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = FALSE, cst.mu = FALSE, fix.eps = TRUE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 4c:", cond))
      list(model = "Birth-death constant extinction fraction (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA, alpha = NA, eps = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 4c:", cond))
    },
    finally = message("Model 4c completed")
  )  

  # Model 4d
  m4d = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = FALSE, cst.mu = FALSE, fix.eps = FALSE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 4d:", cond))
      list(model = "Birth-death varying speciation and extinction (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA, alpha = NA, mu0 = NA, beta = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 4d:", cond))
    },
    finally = message("Model 4d completed")
  )
   
  # Model 5
  m5 = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = TRUE, mu.0 = TRUE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 5:", cond))
      list(model = "Pure birth constant speciation (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 5:", cond))
    },
    finally = message("Model 5 completed")
  )   
  
  # Model 6
  m6 = tryCatch(
    {
      fit_coal_var(tree, lamb0 = 0.01, alpha = -0.001, mu0 = 0.0, beta = 0,
                   cst.lamb = FALSE, mu.0 = TRUE)
    }, 
    error = function(cond) {
      message(paste("Error in Model 6:", cond))
      list(model = "Pure birth varying speciation (coalescent approx)",
           LH = NA, aicc = NA, lamb0 = NA, alpha = NA)
    },
    warning = function(cond) {
      message(paste("Warning in Model 6:", cond))
    },
    finally = message("Model 6 completed")
  )   

  out = data.frame(model = c('1', '2', '3', '4a', '4b', '4c', '4d', '5', '6'),
                   name = c(m1$model, m2$model, m3$model, m4a$model, m4b$model, 
                            m4c$model, m4d$model, m5$model, m6$model),
                   LH = c(m1$LH, m2$LH, m3$LH, m4a$LH, m4b$LH, 
                            m4c$LH, m4d$LH, m5$LH, m6$LH),
                   aicc = c(m1$aicc, m2$aicc, m3$aicc, m4a$aicc, m4b$aicc, 
                            m4c$aicc, m4d$aicc, m5$aicc, m6$aicc),
                   tau0 = c(m1$tau0, m2$tau0, rep(NA, 7)),
                   gamma = c(NA, m2$gamma, rep(NA, 7)),
                   lamb0 = c(NA, NA, m3$lamb0, m4a$lamb0, m4b$lamb0, 
                            m4c$lamb0, m4d$lamb0, m5$lamb0, m6$lamb0),
                   mu0 = c(NA, NA, m3$mu0, m4a$mu0, m4b$mu0, 
                             NA, m4d$mu0, NA, NA),
                   alpha = c(NA, NA, NA, m4a$alpha, NA, m4c$alpha, m4d$alpha, 
                             NA, m6$alpha),
                   beta = c(NA, NA, NA, NA, m4b$beta, NA, m4d$beta, NA, NA),
                   eps = c(rep(NA, 5), m4c$eps, rep(NA, 3)))

  return(out)                 

}

# Fit the 9 models from Morlon et al. 2010 to the 4 diversification scenarios
# from Hurlbert & Stegen 2014, Frontiers in Genetics

#Energy gradient
t4065.30k = read.tree('z:/git/bamm-simulations/sim4065-30k/extant_phy4065_30k.tre')
#Speciation gradient
t5525 = read.tree('z:/git/bamm-simulations/sim5525/extant_phy5525.tre')
t5525.30k = read.tree('z:/git/bamm-simulations/sim5525-30k/run7/extant_phy5525_30k.tre')
#Disturbance gradient
t3865 = read.tree('z:/git/bamm-simulations/sim3865/extant_phy3865.tre')
#Niche conservatism
t3465 = read.tree('z:/git/bamm-simulations/sim3465/run3/extant_phy3465.tre')

# Scale trees to have max branch length of 100
# (advice from Dan Rabosky to get BAMM to fit trees better without getting stuck
# on local optima)
t4065.30ksc = t4065.30k
t4065.30ksc$edge.length = 100*t4065.30ksc$edge.length/max(t4065.30ksc$edge.length)

t5525sc = t5525
t5525sc$edge.length = 100*t5525sc$edge.length/max(t5525sc$edge.length)

t5525.30ksc = t5525.30k
t5525.30ksc$edge.length = 100*t5525.30ksc$edge.length/max(t5525.30ksc$edge.length)

t3865sc = t3865
t3865sc$edge.length = 100*t3865sc$edge.length/max(t3865sc$edge.length)

panda4065.30k = multi.panda.fit(t4065.30ksc)
panda5525sc = multi.panda.fit(t5525sc)
panda3865 = multi.panda.fit(t3865sc)
panda5525.30k = multi.panda.fit(t5525.30ksc)
panda3465 = multi.panda.fit(t3465)

combined = rbind(panda4065.30k, panda5525sc, panda5525.30k, panda3865, panda3465)
combined$sim.id = rep(c('4065-30k', '5525', '5525-30k', '3865', '3465'), each = 9)
combined = combined[, c(12, 1:11)]
combined[,6:12] = signif(combined[, 6:12], 3)

write.csv(combined, 'analysis_output/RPANDA_analysis/panda_output.csv', row.names=F)
