make_sim_matrix = function(curr.matrix,reg.of.origin,w,alpha,alpha.to.beta,sigma_E,carry.cap,energy.gradient,max.K,num.of.bins,max.time,max.richness,replicates) {

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
            for (replicates.use in replicates) {

					curr.matrix = rbind(curr.matrix,c(max(as.numeric(curr.matrix$sim.id)) +1,'to.run',origin,w.use,alpha.use,get.beta*alpha.use,sigma.use,carry.use,gradient.use,max.K.use,bins.use,max.time.use,max.richness.use,replicates.use));
            
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