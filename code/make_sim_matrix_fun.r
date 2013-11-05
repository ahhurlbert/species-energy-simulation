make_sim_matrix = function(curr.matrix,reg.of.origin,w,alpha,alpha.to.beta,sigma_E,carry.cap,energy.gradient,max.K,num.of.bins,max.time,max.richness,disturb_frequency,temperate_disturb_intensity,tropical_disturb_intensity,replicates) {

	for (origin in reg.of.origin) {
	  for (w.use in w) {
	    for (alpha.use in alpha) {
		    for (get.beta in alpha.to.beta) {
		      for (sigma.use in sigma_E) {
		        for (carry.use in carry.cap) {
			        for (gradient.use in energy.gradient) {
			          for (max.K.use in max.K) {
			            for (bins.use in num.of.bins) {
				            for (max.time.use in max.time) {
				              for (max.richness.use in max.richness) {
                        for (disturb_frequency.use in disturb_frequency) {
                          for (temperate_disturb_intensity.use in temperate_disturb_intensity) {
                            for (tropical_distrub_intensity.use in tropical_disturb_intensity) {
                              for (replicates.use in replicates) {

					                      curr.matrix = rbind(curr.matrix,c(max(as.numeric(curr.matrix$sim.id)) +1,'to.run',origin,w.use,alpha.use,get.beta*alpha.use,sigma.use,carry.use,gradient.use,max.K.use,bins.use,max.time.use,max.richness.use,replicates.use,disturb_frequency.use,temperate_disturb_intensity.use,tropical_distrub_intensity.use));
            
                              }
                            }
                          }
                        }
				              }
				            }
			            }
			          }
			        }
		        }
		      }
		    }
	    }
	  }
	}

	return(curr.matrix)
	
};