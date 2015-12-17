# Code for analyzing the simulated phylogenies from 
# Hurlbert & Stegen 2014a,b that were generated under
# different eco-evolutionary dynamics using Morlon et al.'s
# (2010, PLoS Biology) RPANDA package

# Set working directory to species-energy-simulation repo
setwd('z:/git/species-energy-simulation')

source('code/unzipping_files.r')
source('code/time.slice.phylo.r')
library(RPANDA)
library(geiger)
library(paleotree)

# Function for rescaling the branch lengths, with a default setting
# the maximum branch length to 100
rescaleBranchLengths = function(tree, maxLength = 100) {
  tree.out = tree
  tree.out$edge.length = maxLength * tree$edge.length/max(tree$edge.length)
  return(tree.out)
}


# If append = TRUE, then results will be appended to data.frame
# called prevOutput which must have the same column names
multi.panda.fit = function(simID, tree, scale = TRUE, 
                           append = FALSE, prevOutput, write = TRUE) {
  
  message(paste("Fitting PANDA on sim", simID))
  
  if(scale) {
    tree = rescaleBranchLengths(tree)
  }
  
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
      message("Trying again with gamma = 1")
      tryCatch(
        {
          fit_coal_cst(tree, tau0 = 1e-4, gamma = 1, cst.rate = FALSE)
        },
        error = function(cond) {
          message(paste("Error in Model2:", cond))
          list(model = "Equilibrium variable rate", LH = NA, aicc = NA, tau0 = NA, gamma = NA)
        },
        warning = function(cond) {
          message(paste("Warning in Model 2:", cond))
        }
      )
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

  aiccs = c(m1$aicc, m2$aicc, m3$aicc, m4a$aicc, m4b$aicc, m4c$aicc,
            m4d$aicc, m5$aicc, m6$aicc)
  delta.aicc = aiccs - min(aiccs, na.rm = T)

  out = data.frame(sim.id = rep(simID, 9),
                   model = c('1', '2', '3', '4a', '4b', '4c', '4d', '5', '6'),
                   name = c(m1$model, m2$model, m3$model, m4a$model, m4b$model, 
                            m4c$model, m4d$model, m5$model, m6$model),
                   LH = c(m1$LH, m2$LH, m3$LH, m4a$LH, m4b$LH, 
                            m4c$LH, m4d$LH, m5$LH, m6$LH),
                   aicc = c(m1$aicc, m2$aicc, m3$aicc, m4a$aicc, m4b$aicc, 
                            m4c$aicc, m4d$aicc, m5$aicc, m6$aicc),
                   delta.aicc = aiccs - min(aiccs, na.rm = T),
                   w = round(exp(-0.5*delta.aicc)/sum(exp(-0.5*delta.aicc)), 2),
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

  if(append) {
    out = rbind(prevOutput, out)
  }
  
  if(write) {
    write.csv(out, paste('analysis_output/RPANDA_analysis/panda_output_', 
              Sys.Date(), '.csv', sep = ''), row.names = F)
  }
  
  message(paste("Sim", simID, "completed at", Sys.time()))
  
  return(out)                 
}

# This function only works 
getExtantTree = function(simID) {
  simoutput = output.unzip('archived_sim_output', simID)
  phy = simoutput$phylo.out
  all.pops = simoutput$all.populations
  max.time = max(all.pops$time.of.origin)
  extant.pops = subset(all.pops, extant==1)
  extant.phy = drop.tip(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
  output = list(tree = extant.phy, t = max.time, all.pops = all.pops)
  return(output)
}



# Important to set stringsAsFactors to FALSE or else new simID's will prevent
# rbinding
prevOutput = read.csv('z:/git/species-energy-simulation/analysis_output/RPANDA_analysis/panda_output.csv', 
                      header=T, stringsAsFactors = FALSE)

# Fit the 9 models from Morlon et al. 2010 to the 4 diversification scenarios
# from Hurlbert & Stegen 2014, Frontiers in Genetics

#Energy gradient
t4065.30k = read.tree('z:/git/bamm-simulations/sim4065-30k/extant_phy4065_30k.tre')
t4065 = read.tree('z:/git/bamm-simulations/sim4065/run6/extant_phy4065_scaled.tre')
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
t4065.30ksc = rescaleBranchLengths(t4065.30k)
t5525sc = rescaleBranchLengths(t5525)
t5525.30ksc = rescaleBranchLengths(t5525.30k)
t3865sc = rescaleBranchLengths(t3865)

panda4065.30k = multi.panda.fit('4065-30k', t4065.30ksc)
panda5525sc = multi.panda.fit('5525', t5525sc)
panda3865 = multi.panda.fit('3865', t3865sc)
panda5525.30k = multi.panda.fit('5525-30k', t5525.30ksc)
panda3465 = multi.panda.fit('3465', t3465)
panda4065 = multi.panda.fit('4065', t4065)


combined = rbind(panda4065.30k, panda5525sc, panda5525.30k, panda3865, panda3465)
combined[,7:14] = signif(combined[, 7:14], 3)
combined$w = round(combined$w, 2)

write.csv(combined, 'analysis_output/RPANDA_analysis/panda_output.csv', row.names=F)



# Workflow for many sims

# Shell commands for unzipping sim output folders
# for ((i=3565; i<=3574; i++)); do unzip senc.out.$i.zip; done
#sims = c(4065:4084, 3465:3474, 3565:3574, 3866:3874, 3965:3974, 5525:5544, 5625:5644)
sims = 5525:5544

# Read in existing panda output
prevOutput = read.csv('z:/species-energy-simulation/analysis_output/RPANDA_analysis/panda_output.csv', 
                      header=T, stringsAsFactors = FALSE)
# Mid-simulation time point at which to conduct analyses
time = 30000

for (s in sims) {
  # For most sims which are in archived_sim_output folder of github repo
  #simoutput = output.unzip('archived_sim_output', s)
  #phy = simoutput$phylo.out
  #all.pops = simoutput$all.populations
  
  # For sims 5525:5544 which are in particular non-zipped folders
  phy = read.tree(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                        s, '_out/SENC_phylo_sim', s, '.tre', sep = ''))
  all.pops = read.csv(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                            s, '_out/SENC_all.pops_sim', s, '.csv', sep = ''), header=T)
  time.richness = read.csv(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                                 s, '_out/SENC_time.rich_sim', s, '.csv', sep = ''), header=T)
  simoutput = list(all.populations = all.pops, phylo.out = phy, time.richness = time.richness)
  
  max.time = max(all.pops$time.of.origin)
  extant.pops = subset(all.pops, extant==1)
  extant.phy = drop.tip(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
  
  newOutput = multi.panda.fit(simID = s, extant.phy, scale = TRUE, append = TRUE, 
                              prevOutput, write = TRUE)  
  
  if (max.time > time) {
    timeSlice = time.slice.phylo(simoutput, time)
    
    newerOutput = multi.panda.fit(simID = paste(s, "-", time, sep = ""), 
                                  timeSlice$slicedphylo,scale = TRUE, 
                                  append = TRUE, newOutput, write = TRUE)
    prevOutput = newerOutput
  } else {
    prevOutput = newOutput
  }
}

#-----------ANALYZE PANDA OUTPUT------------------------------------------------
panda = read.csv('analysis_output/RPANDA_analysis/panda_output_2015-12-15.csv', header=T)
simkey = read.csv('analysis_output/RPANDA_analysis/simkey.csv', header=T)

panda2 = merge(panda, simkey, by = 'sim.id', all.x = T)

mean.w = aggregate(panda2$w, by = list(panda2$scenario, panda2$origin, panda2$time, panda2$model), 
                   function(x) mean(x, na.rm = T))
names(mean.w) = c('scenario', 'origin', 'time', 'model', 'w')

mean.w = mean.w[order(mean.w$scenario, mean.w$origin, mean.w$time),]

mean.w$col = 'salmon'
mean.w$col[mean.w$scenario == "energy gradient"] = 'limegreen'
mean.w$col[mean.w$scenario == "speciation gradient"] = 'mediumslateblue'
mean.w$col[mean.w$scenario == "pure niche conservatism"] = 'gray50'

scenarios = c('disturbance gradient', 'energy gradient', 'speciation gradient', 
              'pure niche conservatism')

# Plotting mean Akaike weights
pdf('analysis_output/RPANDA_analysis/panda_model_weights.pdf', height = 8, width = 11)
par(mfrow = c(3, 4), mar = c(3,3,3,1), oma = c(4, 4, 0, 0))
for (s in scenarios) {
  for (o in c('temperate', 'tropical')) {
    temp = subset(mean.w, scenario == s & origin == o)
    for (t in unique(temp$time)) {
      temp2 = subset(temp, time == t)
      barplot(temp2$w, ylim = c(0, 1), names.arg = temp2$model, cex.main = 0.9,
              main = paste(s, ", ", o, " origin, \ntime = ", t, sep = ""),
              col = temp2$col[1])
    }
  }
}
mtext("Akaike weight", 2, outer=T, cex = 2, line = 1.75)
mtext("Model", 1, outer = T, cex = 2, line = 2)
dev.off()