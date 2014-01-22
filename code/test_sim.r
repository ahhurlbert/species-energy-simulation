# Use this function to test that the senc_sim_fun function is working properly

test_sim = function() {

  require(ape)
  
    sim_matrix = data.frame(sim.id = 0,
                          status = 'test',
                          reg.of.origin = 'tropical',
                          w = 3,
                          alpha = 0.0001,
                          beta = 0.01, 
                          sigma_E = 1,
                          carry.cap = 'on',
                          energy.gradient = 'on',
                          max.K = 40000,
                          num.of.bins = 11,
                          max.time = 100,
                          max.richness = 3000,
                          replicate = 0,
                          disturb_frequency = 0,
                          temperate_disturb_intensity = 0,
                          tropical_disturb_intensity = 0)
  
  seed = 999
  set.seed(seed)
  sim.results = senc_sim_fun(sim_matrix, sim = 0)
      
  tot.uniq.spp = length(unique(sim.results$all.populations$spp.name))
  extant.spp = nrow(unique(sim.results$all.populations[ sim.results$all.populations$extant == 1, c('spp.name','extant')]))
  time.rich = sim.results$time.richness[sim.results$time.richness$time==100 & 
                                          !sim.results$time.richness$region %in% c(0,11),]
  gamma = gammaStat(sim.results$phylo)
  
  if(tot.uniq.spp == 1945) { tot.uniq.check = 'PASS' } else { tot.uniq.check = 'FAIL' }
  if(extant.spp == 1799) { extant.spp.check = 'PASS' } else { extant.spp.check = 'FAIL' }
  if(ceiling(gamma*1e5)/1e5 == -22.27503) { gamma.check = 'PASS' } else { gamma.check = 'FAIL' }
  if( identical(sim.results$time.richness$spp.rich[1190:1199], c(2, 47, 106, 174, 293, 359, 473, 565, 663, 713)) ) {
    time.rich.check = 'PASS'
  } else { time.rich.check = 'FAIL' }
  
  cat("Simulation test results:", paste("  Total unique species:", tot.uniq.check),
      paste("  Total extant species:", extant.spp.check),
      paste("  Richness gradient at t = 100:", time.rich.check),
      paste("  Gamma:", gamma.check), sep = "\n")
  
} #end function

