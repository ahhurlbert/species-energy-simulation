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
prevOutput = read.csv('z:/git/species-energy-simulation/analysis_output/RPANDA_analysis/panda_output_2015-12-17.csv', 
                      header=T, stringsAsFactors = FALSE)

# Fit the 9 models from Morlon et al. 2010 to the 4 diversification scenarios
# from Hurlbert & Stegen 2014, Frontiers in Genetics

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
  if (!s %in% c(5525:5544, 5625:5644)) {
    simoutput = output.unzip('archived_sim_output', s)
    phy = simoutput$phylo.out
    all.pops = simoutput$all.populations
  }
  # For sims from Frontiers paper
  if (s %in% c(5525:5544, 5625:5644)) {
    phy = read.tree(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                          s, '_out/SENC_phylo_sim', s, '.tre', sep = ''))
    all.pops = read.csv(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                              s, '_out/SENC_all.pops_sim', s, '.csv', sep = ''), header=T)
    time.richness = read.csv(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                                   s, '_out/SENC_time.rich_sim', s, '.csv', sep = ''), header=T)
    simoutput = list(all.populations = all.pops, phylo.out = phy, time.richness = time.richness)
                          
  }
    
  max.time = max(all.pops$time.of.origin)
  extant.pops = subset(all.pops, extant==1)
  extant.phy = drop.tip(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
  
  newOutput = prevOutput
  #newOutput = multi.panda.fit(simID = s, extant.phy, scale = TRUE, append = TRUE, 
  #                            prevOutput, write = TRUE)  
  
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

#######################################

# Fit the Manceau et al. 2015 SGD model

# (which in the future should be added to multi.panda.fit)

# Currently doesn't work (at least on sims 4065, 4066), gives this error:

# DLSODA-  TOUT(=R1) too close to T(=R2) to start integration.
# In above message, R1 = 100.281, R2 = 100.281

#Error in lsoda(y, times, func, parms, ...) : 
#  illegal input detected before taking any integration steps - see written message 


# Mid-simulation time point at which to conduct analyses
time = 30000

for (s in sims) {
  # For most sims which are in archived_sim_output folder of github repo
  if (!s %in% c(5525:5544, 5625:5644)) {
    simoutput = output.unzip('archived_sim_output', s)
    phy = simoutput$phylo.out
    all.pops = simoutput$all.populations
  }
  # For sims from Frontiers paper
  if (s %in% c(5525:5544, 5625:5644)) {
    phy = read.tree(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                          s, '_out/SENC_phylo_sim', s, '.tre', sep = ''))
    all.pops = read.csv(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                              s, '_out/SENC_all.pops_sim', s, '.csv', sep = ''), header=T)
    time.richness = read.csv(paste('z:/manuscripts/frontierstropicaldiversity/raw_sim_output/sim',
                                   s, '_out/SENC_time.rich_sim', s, '.csv', sep = ''), header=T)
    simoutput = list(all.populations = all.pops, phylo.out = phy, time.richness = time.richness)
    
  }
  
  max.time = max(all.pops$time.of.origin)
  extant.pops = subset(all.pops, extant==1)
  extant.phy = drop.tip(phy, phy$tip.label[!phy$tip.label %in% extant.pops$spp.name])
  tree = rescaleBranchLengths(extant.phy)
  
  tot_time = max(node.age(tree)$ages)
  par_init = c(1e7, 1e7-0.5, 1)  #initial parameters as used in documentation example
  sgd_out = fit_sgd(tree, tot_time, par_init, f = 1) 
  
  out = data.frame(sim.id = s,
                   model = 'sgd',
                   name = 'Speciation by Genetic Differentiation',
                   LH = sgd_out$LH,
                   aicc = sgd_out$aicc,
                   delta.aicc = NA,
                   w = NA,
                   tau0 = NA,
                   gamma = NA,
                   lamb0 = NA,
                   mu0 = NA,
                   alpha = NA,
                   beta = NA,
                   eps = NA, 
                   birth = sgd_out$par$birth,
                   growth = sgd_out$par$growth,
                   mutation = sgd_out$par$mutation)
  
  prevOutput = rbind(prevOutput, out)
  if (max.time > time) {
    timeSlice = time.slice.phylo(simoutput, time)
    tree2 = rescaleBranchLengths(timeSlice$slicedphylo)
    tot_time2 = max(node.age(tree2)$ages)
    par_init2 = c(1e7, 1e7-0.5, 1)
    sgd_out2 = fit_sgd(tree2, tot_time2, par_init2, f = 1) 
    
    out2 = data.frame(sim.id = paste(s, '-', time, sep = ''),
                     model = 'sgd',
                     name = 'Speciation by Genetic Differentiation',
                     LH = sgd_out2$LH,
                     aicc = sgd_out2$aicc,
                     delta.aicc = NA,
                     w = NA,
                     tau0 = NA,
                     gamma = NA,
                     lamb0 = NA,
                     mu0 = NA,
                     alpha = NA,
                     beta = NA,
                     eps = NA, 
                     birth = sgd_out2$par$birth,
                     growth = sgd_out2$par$growth,
                     mutation = sgd_out2$par$mutation)
    
    prevOutput = rbind(prevOutput, out2)
  }
  message(paste("Sim", s, "completed at", Sys.time()))
}



#-----------ANALYZE PANDA OUTPUT------------------------------------------------
panda = read.csv('analysis_output/RPANDA_analysis/panda_output_2015-12-17.csv', header=T)
simkey = read.csv('analysis_output/RPANDA_analysis/simkey.csv', header=T)
modelkey = read.csv('analysis_output/RPANDA_analysis/modelkey.csv', header=T)

panda2 = merge(panda, simkey, by = 'sim.id', all.x = T)
panda3 = merge(panda2, modelkey, by = 'model', all.x = T)

# Mean, SD, and SE of Akaike weights by original Morlon et al. 2010 model number
w.summary = aggregate(panda3$w, by = list(panda3$scenario, panda3$origin, panda3$time, 
                                          panda3$model, panda3$modelgrp), 
                      function(x) c(mean = mean(x, na.rm = T),
                                     sd = sd(x, na.rm = T),
                                     n = length(x[!is.na(x)])))
w.summary = do.call(data.frame, w.summary)
names(w.summary) = c('scenario', 'origin', 'time', 'model', 'modelgrp', 'w.mean', 'w.sd', 'w.n')
w.summary$w.se = w.summary$w.sd/sqrt(w.summary$w.n)

w.summary = w.summary[order(w.summary$scenario, w.summary$origin, w.summary$time),]

# Mean, SD, and SE of Akaike weights by model group
# --a: saturated diversity (Models 1-2)
# --b: expanding diversity with positive extinction (Models 3-4)
# --c: expanding diversity with no extinction (Models 5+6)
w.groups = aggregate(panda3$w, by = list(panda3$sim.id, panda3$scenario, panda3$origin, panda3$modelgrp, panda3$time), 
                      function(x) sum(x, na.rm = T))
names(w.groups) = c('sim.id', 'scenario', 'origin', 'modelgrp', 'time', 'w')
w.groupsumm = aggregate(w.groups$w, by = list(w.groups$scenario, w.groups$origin, w.groups$time, w.groups$modelgrp), 
                        function(x) c(mean = mean(x, na.rm = T),
                                      sd = sd(x, na.rm = T),
                                      n = length(x[!is.na(x)])))
w.groupsumm = do.call(data.frame, w.groupsumm)
names(w.groupsumm) = c('scenario', 'origin', 'time', 'model', 'w.mean', 'w.sd', 'w.n')
w.groupsumm$w.se = w.groupsumm$w.sd/sqrt(w.groupsumm$w.n)

w.groupsumm = w.groupsumm[order(w.groupsumm$scenario, w.groupsumm$origin, w.groupsumm$time),]





w.summary$col = 'salmon'
w.summary$col[w.summary$scenario == "energy gradient"] = 'limegreen'
w.summary$col[w.summary$scenario == "speciation gradient"] = 'mediumslateblue'
w.summary$col[w.summary$scenario == "pure niche conservatism"] = 'gray50'

w.summary$col2 = 'salmon'
w.summary$col2[w.summary$modelgrp == "b"] = 'limegreen'
w.summary$col2[w.summary$modelgrp == "c"] = 'mediumslateblue'




scenarios = c('energy gradient', 'speciation gradient', 'disturbance gradient',
              'pure niche conservatism')

# Plotting mean Akaike weights
pdf('analysis_output/RPANDA_analysis/panda_model_weights.pdf', height = 10, width = 11)
par(mfrow = c(4, 4), mar = c(3,3,3,1), oma = c(4, 4, 0, 0))
for (s in scenarios) {
  for (o in c('temperate', 'tropical')) {
    temp = subset(w.summary, scenario == s & origin == o)
    for (t in unique(temp$time)) {
      temp2 = subset(temp, time == t)
      barCenters = barplot(temp2$w.mean, ylim = c(0, 1.15), names.arg = temp2$model, 
              main = paste(s, ", ", o, " origin, \ntime = ", t, sep = ""),
              cex.main = 0.9, col = temp2$col[1], las = 1)
      arrows(barCenters, temp2$w.mean - temp2$w.se*2, barCenters, 
             temp2$w.mean + temp2$w.se*2, angle = 90, code = 3, length = 0.03)
    }
  }
}
mtext("Akaike weight", 2, outer=T, cex = 2, line = 1.75)
mtext("Model", 1, outer = T, cex = 2, line = 2)
dev.off()

# Plot of tropical origin sims at 2 time points
salmonbg = c(254, 230, 226)
greenbg = c(214, 245, 214)
bluebg = c(229, 225, 252)

pdf('analysis_output/RPANDA_analysis/panda_model_weights_tropical.pdf', 
    height = 6, width = 11)
par(mfcol = c(2, 4), mar = c(3,3,2,1), oma = c(4, 4, 4, 0))
for (s in scenarios) {
  temp = subset(w.summary, scenario == s & origin == 'tropical')
  for (t in unique(temp$time)) {
    temp2 = subset(temp, time == t)
    barCenters = barplot(rep(0,9), ylim = c(0, 1.15), las = 1, cex.axis = 1.25)
    width = barCenters[2] - barCenters[1]
    rect(0,0, barCenters[2] + width/2, 1.1, border = NA, 
         col = rgb(salmonbg[1], salmonbg[2], salmonbg[3], maxColorValue = 255))
    rect(barCenters[2] + width/2, 0, barCenters[7] + width/2, 1.1, border = NA,
         col = rgb(greenbg[1], greenbg[2], greenbg[3], maxColorValue = 255))
    rect(barCenters[7] + width/2, 0, barCenters[9] + width/2, 1.1, border = NA,
         col = rgb(bluebg[1], bluebg[2], bluebg[3], maxColorValue = 255))
    barCenters = barplot(temp2$w.mean, ylim = c(0, 1.15), xaxt = "n", yaxt = "n",
                         col = temp2$col2, add = T)
    arrows(barCenters, temp2$w.mean - temp2$w.se*2, barCenters, 
           temp2$w.mean + temp2$w.se*2, angle = 90, code = 3, length = 0.03)
    axis(1, temp2$model[seq(1,9,2)], at = barCenters[seq(1,9,2)], cex.axis = 1.2)
    axis(1, temp2$model[seq(2,9,2)], at = barCenters[seq(2,9,2)], cex.axis = 1.2)
  }
}
mtext("Akaike weight", 2, outer=T, cex = 2, line = 1.75)
mtext("Model", 1, outer = T, cex = 2, line = 2)
mtext(c("Energy\ngradient", "Speciation\ngradient", "Disturbance\ngradient", 
        "Pure niche\nconservatism"), 3, at = c(1/8, 3/8, 5/8, 7/8), outer = T, cex = 1.75, 
      line = -1)
dev.off()

