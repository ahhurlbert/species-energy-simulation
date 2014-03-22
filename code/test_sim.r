# Use this function to test that the senc_sim_fun function is working properly

test_sim = function() {

  require(ape)
  source('code/senc_sim_fun.r')
  source('code/supplemental_analysis_functions.r')
    
  sim_matrix = data.frame(sim.id = 99999,
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
                          tropical_disturb_intensity = 0,
                          specn.gradient = 'off',
                          specn.factor = NA,
                          gamma = 0.1)
  
  sim.results = senc_sim_fun(sim_matrix, sim = 99999)
      
  tot.uniq.spp = length(unique(sim.results$all.populations$spp.name))
  extant.spp = nrow(unique(sim.results$all.populations[ sim.results$all.populations$extant == 1, c('spp.name','extant')]))
  time.rich = sim.results$time.richness[sim.results$time.richness$time==100 & 
                                          !sim.results$time.richness$region %in% c(0,11),]
  gamma = gammaStat(sim.results$phylo)
  
  if(tot.uniq.spp == 1647) { tot.uniq.check = 'PASS' } else { tot.uniq.check = 'FAIL' }
  if(extant.spp == 1449) { extant.spp.check = 'PASS' } else { extant.spp.check = 'FAIL' }
  if(ceiling(gamma*1e5)/1e5 == -18.42047) { gamma.check = 'PASS' } else { gamma.check = 'FAIL' }
  if( identical(sim.results$time.richness$spp.rich[1190:1199], c(0, 0, 6, 51, 85, 191, 370, 513, 693, 760)) ) {
    time.rich.check = 'PASS'
  } else { time.rich.check = 'FAIL' }
  
  cat("Simulation test results:", paste("  Total unique species:", tot.uniq.check),
      paste("  Total extant species:", extant.spp.check),
      paste("  Richness gradient at t = 100:", time.rich.check),
      paste("  Gamma:", gamma.check), sep = "\n")
  
} #end function

