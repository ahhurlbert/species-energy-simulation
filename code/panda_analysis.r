# Code for analyzing the simulated phylogenies from 
# Hurlbert & Stegen 2014a,b that were generated under
# different eco-evolutionary dynamics using Morlon et al.'s
# (2010, PLoS Biology) RPANDA package

# Set working directory to species-energy-simulation repo

library(RPANDA)
library(geiger)

# Read in trees
sim1 = 4065
sim2 = 3465

tree1 = read.tree('raw_sim_output/sim3465_out/SENC_phylo_sim4065.tre')
tree2 = read.tree('raw_sim_output/sim3465_out/SENC_phylo_sim3465.tre')

# Prune out extinct taxa 
# Says there are 2610 species, but original sim analysis output says should be 2212
# --I think this difference is because in the sim analysis, species that occur only
#   in the boundary regions (0 or 11) get dropped. Shouldn't matter for this
#   analysis whether they are included or not, I wouldn't think.
tree1p = drop.extinct(tree1) #took ~2 hrs
tree2p = drop.extinct(tree2) #took ~15 min

# Using sim analysis code:
sim_dir = 'z:/SENCoutput/Hurlbert_and_Stegen_2014/raw_sim_output/senc.out.4065'
all.populations = read.csv(paste(sim_dir,'/SENC_all.pops_sim',sim,'.csv',sep=''), header=T)
time.richness = read.csv(paste(sim_dir,'/SENC_time.rich_sim',sim,'.csv',sep=''), header=T)
phylo.out = read.tree(paste(sim_dir,'/SENC_phylo_sim',sim,'.tre',sep=''))
params.out = read.csv(paste(sim_dir,'/SENC_params.out_sim',sim,'.csv',sep=''), header=T)

t=params.out$max.time
sub.species = as.character(unique(subset(all.populations,time.of.sp.origin <= t & time.of.sp.extinction > t, select = 'spp.name'))[,1]);
# Some species may be extant globally (extant==1) but in our boundary regions (0,11) only;
# we need to eliminate species that are not extant within regions 1-10 (which is all that is
# reflected in the all.populations dataframe)
time.slice.populations = all.populations;
time.slice.populations$extant = 0;
time.slice.populations$extant[time.slice.populations$time.of.origin <= t & time.slice.populations$time.of.extinction > t] = 1
extant.ornot = aggregate(time.slice.populations$extant,by=list(time.slice.populations$spp.name),sum)
extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
sub.species2 = sub.species[!sub.species %in% extinct.species]
tips.to.drop = as.character(phylo.out$tip.label[!phylo.out$tip.label %in% sub.species2]);
sub.phylo = drop.tip(phylo.out,tips.to.drop);
sub.phylo$Nnode



# Following doesn't work with tau0 = .01, try other values besides 1e-4
result1p = fit_coal_cst(tree1.2, tau0 = 1e-4, gamma = 1, cst.rate = TRUE)

# This takes ~2 hours:
result2 <- fit_coal_var(tree2, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0)
# This one is instantaneous:
result2.cst = fit_coal_cst(tree2, tau0 = 0.01, gamma = 1, cst.rate = TRUE)

# This one results in an error:
# Error in optim(init, optimLH.coalBD, method = meth, control = list(ndeps = 10^(-4))) : 
# function cannot be evaluated at initial parameters
result1 <- fit_coal_var(tree1, lamb0=0.01, alpha=-0.001, mu0=0.0, beta=0)


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




write.csv(out, 'analysis_output/RPANDA_analysis/sim3465_panda.csv', row.names=F)































# Per RPANDA documentation, 

# Use fit_coal_cst() for Models 1 and 2 and fit_coal_var for Models 3-6 from Morlon et al. 2010


output = c()
for (i in sims) {
  treename = paste('raw_sim_output/sim', i, '_out/SENC_phylo_sim', i, '.tre', sep = '')
  tree = read.tree(treename)
  out = fit_coal_cst(tree, tau0=1.e-2, gamma=-1, cst.rate=TRUE, N0 = length(tree$tip.label))
}


fit_coal_cst = function (phylo, tau0 = 0.01, gamma = 1, cst.rate = FALSE, meth = "Nelder-Mead", 
          N0 = 0) 
{
  if (!inherits(phylo, "phylo")) 
    stop("object \"phylo\" is not of class \"phylo\"")
  Vtimes <- sort(branching.times(phylo))
  ntips <- Ntip(phylo)
  if (N0 == 0) {
    N0 <- ntips
  }
  if (cst.rate == FALSE) {
    init <- c(tau0, gamma)
    nbpar <- length(init)
    nbobs <- length(Vtimes) - 1
    optimLH <- function(init) {
      tau0 <- init[1]
      gamma <- init[2]
      LH <- likelihood_coal_cst(Vtimes, ntips, tau0, gamma, 
                                N0)$res
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    res <- list(model = "Equilibrium variable rate", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * nbpar + 2 * nbpar * (nbpar + 
                                                                   1)/(nbobs - nbpar - 1), tau0 = temp$par[1], gamma = temp$par[2])
    return(res)
  }
  else {
    init <- c(tau0)
    nbpar <- length(init)
    nbobs <- length(Vtimes) - 1
    optimLH <- function(init) {
      tau0 <- init[1]
      LH <- likelihood_coal_cst(Vtimes, ntips, tau0, 0, 
                                N0)$res
      return(-LH)
    }
    temp <- suppressWarnings(optim(init, optimLH, method = meth))
    res <- list(model = "Equilibrium constant rate", LH = -temp$value, 
                aicc = 2 * temp$value + 2 * nbpar + 2 * nbpar * (nbpar + 
                                                                   1)/(nbobs - nbpar - 1), tau0 = temp$par[1])
    return(res)
  }
}

#####
> likelihood_coal_cst
function (Vtimes, ntips, tau0, gamma, N0) 
{
  if (gamma == 0) {
    return(.likelihood_coal_cst_mod(Vtimes, ntips, tau0, 
                                    N0))
  }
  else {
    return(.likelihood_coal_exp_mod(Vtimes, ntips, tau0, 
                                    gamma, N0))
  }
}