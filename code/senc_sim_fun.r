## Simulation function

senc_sim_fun = function(sim.matrix, sim) {

  set.seed(sim)
  
  # Assign simulation parameters from the sim.matrix based on the sim.id
  # --------------------------------------------------------------------
  region.of.origin = sim.matrix$reg.of.origin[sim.matrix$sim.id == sim]  #'temperate' or 'tropical'
  
  w = sim.matrix$w[sim.matrix$sim.id == sim]  #strength of environmental filtering (smaller values = stronger)
  
  alpha = sim.matrix$alpha[sim.matrix$sim.id == sim]  #per capita speciation rate
  
  beta = sim.matrix$beta[sim.matrix$sim.id == sim]  #per capita dispersal rate to adjacent region
  
  gamma = sim.matrix$gamma[sim.matrix$sim.id == sim]  #exponential rate of decay in extinction probability with population size
  
  sigma_E = sim.matrix$sigma_E[sim.matrix$sim.id == sim]  #strength of niche conservatism (smaller values = stronger)
  
  carry.cap = sim.matrix$carry.cap[sim.matrix$sim.id == sim]  #'on' if limits exist to community abundance; else 'off'
  
  energy.gradient = sim.matrix$energy.gradient[sim.matrix$sim.id == sim]  #'on' if maximum abundance varies across the gradient; else 'off'
  
  specn.gradient = sim.matrix$specn.gradient[sim.matrix$sim.id == sim] #'on' if specn rate is higher in tropics; 'off' if does not vary across gradient
  
  specn.factor = sim.matrix$specn.factor[sim.matrix$sim.id == sim] #factor by which tropical speciation rate is greater than temperate
  
  max.K = sim.matrix$max.K[sim.matrix$sim.id == sim]  #num of individuals that can be supported in region with the highest carrying capacity
  
  min.K = max.K/10  #num of individuals that can be supported in region with the lowest carrying capacity when there's a gradient
  
  num.of.bins = sim.matrix$num.of.bins[sim.matrix$sim.id == sim]  #number of spatial bins
  
  max.time = sim.matrix$max.time[sim.matrix$sim.id == sim]  #maximum number of time steps to run simulation
  
  max.richness = sim.matrix$max.richness[sim.matrix$sim.id == sim]  #maximum number of species before simulation breaks off
  
  temperate_disturb_intensity = sim.matrix$temperate_disturb_intensity[sim.matrix$sim.id == sim]
    #fraction of individuals killed in disturbance event at temperate end of gradient
  
  tropical_disturb_intensity = sim.matrix$tropical_disturb_intensity[sim.matrix$sim.id == sim]
    #fraction of individuals killed in disturbance event at tropical end of gradient
  
  disturb_frequency = sim.matrix$disturb_frequency[sim.matrix$sim.id == sim]  #frequency of disturbance in number of time steps
  
  if (disturb_frequency == 0) {
    disturb_times = 0
    reg.disturb.intensity = data.frame(region = 0:num.of.bins, intensity = rep(0, length = num.of.bins + 1)) 
  } else {
    disturb_times = seq(0, max.time, disturb_frequency)
    reg.disturb.intensity = data.frame(region = 0:num.of.bins, 
                                       intensity = seq(temperate_disturb_intensity, tropical_disturb_intensity, length = num.of.bins + 1)) 
  }
  
	if (region.of.origin == 'tropical') {reg.of.origin = num.of.bins - 1}
	if (region.of.origin == 'temperate') {reg.of.origin = 1} 
  #--------------------------------------------------------------------------------------
  
  
  ## Additional fixed and region-specific parameters 
  #-------------------------------------------------
  # abiotic environmental gradient, ranging from 0 to 40 C
  min.env = 0
  max.env = 40
  
  total.mutations = 0 # the cumulative number of mutations

  # regional environments and carrying capacities with and without an energy gradient
  reg.E.K. = data.frame(region = 0:num.of.bins, reg.env = seq(min.env, max.env, length = num.of.bins + 1))
  if (energy.gradient == 'on') { 
    reg.E.K.$carry.cap = seq(min.K, max.K, length = num.of.bins + 1)
  } 
  if (energy.gradient == 'off') { 
    reg.E.K.$carry.cap = rep(max.K, num.of.bins + 1)
  } 
  # regional speciation rates depending on whether they vary across the gradient
  if (specn.gradient == 'on') {
    reg.E.K.$reg.alpha = seq(alpha/(specn.factor ^ .5), alpha*(specn.factor ^ .5), length = num.of.bins + 1)
  }
  if (specn.gradient == 'off') {
    reg.E.K.$reg.alpha = rep(alpha, num.of.bins + 1)
  }
  
  
  ## end additional parameters --------------------------------------------------------
  
  
  ## Initialize output variables, matrices, dataframes
  #------------------------------------------------------------------------------------
	# Set up initial matrix that hold extant populations including their region, name, environmental optimum, time of orgin, 
  # time of extinction, and population size. Also included are the region-level carrying capacity, environmental condition, 
  # and current species richness one matrix for all regions
	all.populations = matrix(-999, nrow = 0, ncol = 8)
	colnames(all.populations) = c('region',
                                'spp.name',
                                'extant',
                                'env.opt',
                                'time.of.origin',
                                'time.of.extinction',
                                'time.of.sp.origin',
                                'time.of.sp.extinction')
  # environmental optimum for ancestor species and initial trait value in the ancestral environment
  all.populations = rbind(all.populations, 
                          c(reg.of.origin, 1, 1, reg.E.K.$reg.env[reg.E.K.$region == reg.of.origin], 0, max.time + 1, 0, max.time + 1)) 
  all.populations = merge(all.populations, reg.E.K., by = 'region')
    
	# set up initial attributes for building the phylogeny
	edge = as.data.frame(rbind(c(1, 2, 1, 1), c(1, 3, 1, 2))) 
  names(edge) = c('from.node', 'to.node', 'alive', 'spp')
	edge.length = rep(NA, 2)
	stem.depth = numeric(2)
	next.node = 4
	next.spp = 2

	# set up initial matrix that will hold time and number of species
	time.richness = matrix(-999, nrow = max.time*(num.of.bins + 1), ncol = 4)
	time.rich.row.id = 1

  #print.times = seq(0,max.time,10)       #used to print optional updates of simulation progress to the console
  
  # set up variables that will tract extinct species/populations
  extinct.pops.output.times = numeric()
  tot.extinct.pops = 0
  # end of output initialization -------------------------------------------------------------------
  
  
  #Simulation process functions (trait mutation, extinction, dispersal, species mutation/speciation)
  #------------------------------------------------------------------------------------------------
	# environmental optimum heritability function, which varies with sigma_E
	mut.opt.fun = function(sigma_E.param, num.to.generate) { 
    rnorm(num.to.generate, mean = 0, sd = sigma_E.param) 
	}
  
  # extinction function a negative exponential function of population size
  extinction.fun = function (gamma.param, population.size) { 
    exp(-gamma.param * population.size) 
	}
	
  # dispersal function
  dispersal.fun = function (population.size, beta.param) { 
    apply(as.matrix(population.size), 1, function(p) sum(rbinom(p, 1, beta.param))) 
	}
  
	# mutation function leading to differentiated species
  mutation.fun = function (population.size, alpha.param) { 
    df = cbind(population.size, alpha.param)
    apply(df, 1, function(x) sum(rbinom(x[1], 1, x[2]))) 
	}

  # environmental fit function used in determining population sizes
  env.fit.fun = function (w.param,E,spp.opt) { 
    exp(-(spp.opt - E)^2/(w.param^2)) 
  }
  
  # function for making a phylogeny (of class 'phylo') from the simulation output
  make.phylo.jimmy.fun = function(t, edge.length.out, edge.out, stem.depth.out) {
    
    edge.length.out[which(edge.out$alive == 1)] = t - stem.depth.out[which(edge.out$alive == 1)]
    edge.out = as.matrix(edge.out)
    n = -1
    for (i in 1:max(edge.out[, c('from.node', 'to.node')])) {
      if (any(edge.out[, 'from.node'] == i)) {
        edge.out[which(edge.out[, 'from.node'] == i), 'from.node'] = n
        edge.out[which(edge.out[, 'to.node'] == i), 'to.node'] = n
        n = n - 1
      }
    }
    
    edge.out[which(edge.out[, c('from.node', 'to.node')] > 0)] = 1:length(edge.out[which(edge.out[, c('from.node', 'to.node')] > 0)])
    tip.label = edge.out[, 'spp'][edge.out[, 'spp'] != -999]
    edge.only = edge.out[, c('from.node', 'to.node')]
    mode(edge.only) = "character"
    mode(tip.label) = "character"
    obj = list(edge = edge.only, edge.length = edge.length.out, tip.label = tip.label)
    class(obj) = "phylo"
    phylo.out = old2new.phylo(obj)
    phylo.out = read.tree(text = write.tree(phylo.out))
    return(phylo.out)
  }
  
  # end simulation functions -----------------------------------------------------------------


	# Start main time loop
  #------------------------------------------------------
	for (curr.time in 1:max.time) {

		## update population sizes
		all.populations$env.fit = env.fit.fun(w.param = w,E = all.populations$reg.env, spp.opt = all.populations$env.opt)
		all.populations$pop.size = all.populations$env.fit * all.populations$carry.cap * all.populations$extant
		reg.pop.sizes = as.data.frame(tapply(all.populations$pop.size, all.populations$region, FUN = sum))
		reg.pop.sizes = cbind(rownames(reg.pop.sizes), reg.pop.sizes) 
    colnames(reg.pop.sizes) = c("region", "reg.individuals")
		all.populations = merge(all.populations, reg.pop.sizes, by = "region")
		
    #The pop.multiplier assures that the sum of species abundances will not exceed carrying capacity
    all.populations$pop.multiplier = all.populations$carry.cap / all.populations$reg.individuals
		all.populations$pop.multiplier[all.populations$pop.multiplier > 1] = 1
    
		# if there is no carrying capacity, then population size is not influenced by the actual 
    # community-level number of individuals, and so the sum across species is not modified (pop.multiplier = 1)
    if (carry.cap == 'off') {all.populations$pop.multiplier = 1} 
		
    all.populations$pop.size = round( all.populations$pop.size * all.populations$pop.multiplier , digits = 0)
		reg.pop.sizes = as.data.frame(tapply(all.populations$pop.size, all.populations$region, FUN = sum))
		reg.pop.sizes = cbind(rownames(reg.pop.sizes), reg.pop.sizes) 
    colnames(reg.pop.sizes) = c("region", "reg.individuals.real")
		all.populations = merge(all.populations, reg.pop.sizes, by = 'region')
    
    # Perturb communities if disturbance has been specified and it's 'time'
    if (is.element(curr.time, disturb_times) == T) {
      for (reg.sampled in unique(all.populations$region[all.populations$extant==1])) {
        reg.disturb = reg.disturb.intensity$intensity[reg.disturb.intensity$region==reg.sampled]
        all.populations$pop.size[all.populations$region == reg.sampled & all.populations$extant==1] =
          round((1 - reg.disturb) * all.populations$pop.size[all.populations$region == reg.sampled & 
                                                                 all.populations$extant == 1], digits = 0)
      }  
    }

		# evaluate mutation across all extant species
		num.of.mutations = mutation.fun(population.size = all.populations$pop.size, alpha.param = all.populations$reg.alpha)
		all.populations$mutate = num.of.mutations

		# find new trait values for mutants
		all.populations$mutant.env.opt = all.populations$env.opt + mut.opt.fun(sigma_E, num.to.generate = nrow(all.populations))

		# evaluate dispersal across all extant species
		all.populations$disperse = dispersal.fun(population.size = all.populations$pop.size, beta.param = beta)
		
		# evaluate extinction across all extant species
		extinction.probs = extinction.fun(gamma.param = gamma, population.size = all.populations$pop.size)
		all.populations$extinction = apply(as.matrix(extinction.probs), 1, function(p) rbinom(1, 1, prob = p))

		# collecting richness data, including only those species that have persisted across at least one time step, 
    # so not including new arrival via dispersal or mutation
		for (i in 0:num.of.bins) {
			time.richness[time.rich.row.id,] = c(curr.time, i, length(all.populations$region[all.populations$region == i & all.populations$extinction==0]),
                                           reg.E.K.$reg.env[reg.E.K.$region==i])
			time.rich.row.id = time.rich.row.id + 1
		}

		# build attributes used in phylogeny construction
		mutated.populations = all.populations[all.populations$mutate > 0,]
		# leave only one mutant from each species. Assuming there are at most two mutations per species, 
    # choose the region of mutation randomly...see the 'sample' piece of the code
		mutated.populations = mutated.populations[is.element(rownames(mutated.populations), 
                                                         rownames(unique(mutated.populations[,c('spp.name','env.opt')],
                                                                         fromLast = sample(c(FALSE,TRUE),1))))==TRUE,]
		if (nrow(mutated.populations) > 0) {
      total.mutations = total.mutations + 1
		}
    # this is just for the first mutation event
		if (nrow(mutated.populations) > 0 & total.mutations == 1) { 
      
			stem.depth = as.numeric(c(curr.time, curr.time))
			all.populations = rbind(all.populations, 
                              c(mutated.populations$region, next.spp, 1, mutated.populations$mutant.env.opt,
                                curr.time, max.time + 1, curr.time, max.time + 1, rep(-999, ncol(all.populations) - 8)))
		}
		
    # this handles all mutations after the first mutation
		if (nrow(mutated.populations) > 0 & total.mutations > 1) { 
		  
			for (i in 1:nrow(mutated.populations)) {

				next.spp = next.spp + 1
				mutated.sp = mutated.populations$spp.name[i] # this is the true species name, not a node name, so it doesn't change.
				mutated.node = edge$to.node[edge$spp == mutated.sp]
				e = edge[which(edge$alive == 1), ]
				edge$alive[which(edge$to.node == mutated.node)] = 0 
				edge$spp[edge$spp == mutated.sp] = -999
				edge = rbind(edge, c(mutated.node, next.node, 1, mutated.sp), c(mutated.node, next.node + 1, 1, next.spp))
				next.node = next.node + 2

				stem.depth <- c(stem.depth, curr.time, curr.time)
				x <- which(edge[, 2] == mutated.node)
				edge.length[x] <- curr.time - stem.depth[x]
				edge.length <- c(edge.length, NA, NA)

				all.populations = rbind(all.populations,
                                c(mutated.populations$region[i], next.spp, 1, mutated.populations$mutant.env.opt[i],
                                  curr.time, max.time + 1, curr.time, max.time + 1, rep(-999, ncol(all.populations)-8)))
			}
		} 

		# disperse species that were selected to do so	
		disp.spp = all.populations[all.populations$disperse > 0, ]
		if ( nrow(disp.spp) > 0 ) {
			disp.spp$region = disp.spp$region + sample(c(-1,1), length(disp.spp$region), replace=T, prob = c(0.5,0.5))
			disp.spp$time.of.origin = curr.time
			disp.spp = disp.spp[is.element(apply(disp.spp[, 1:4], 1, paste, collapse=","),
                                     apply(all.populations[, 1:4], 1, paste, collapse=",")) == F, ]
			all.populations = rbind(all.populations, disp.spp)
			# take out any populations that dispersed beyond the bins being evaluated
			all.populations = all.populations[is.element(all.populations$region, 0:num.of.bins) == TRUE, ] 
		}

		# remove those species that are going extinct
		all.populations$time.of.extinction[all.populations$extinction == 1 & all.populations$extant == 1] = curr.time 
		# --create a vector of species with populations going extinct this time step
		extinct.this.time = unique(all.populations$spp.name[all.populations$extinction == 1 & all.populations$extant == 1]) 
		all.populations$extant[all.populations$extinction == 1 & all.populations$extant == 1] = 0

		# find those species that went globally extinct this time step
		extant.spp = unique(all.populations$spp.name[all.populations$extant == 1])
		global.extinct.this.step = extinct.this.time[is.element(extinct.this.time, extant.spp) == F]
		all.populations$time.of.sp.extinction[is.element(all.populations$spp.name, global.extinct.this.step) == T] = curr.time

		# change the 'alive' status for those species that went globally extinct this time step
		edge$alive[is.element(edge$spp, global.extinct.this.step) == T] = 0
		edge.length[is.element(edge$spp,global.extinct.this.step) == T] = curr.time - 
                                                            stem.depth[is.element(edge$spp, global.extinct.this.step) == T]

		# update columns for carry.cap, reg.env, and reg.richness
		all.populations = all.populations[ , 1:8]
		all.populations = merge(all.populations, reg.E.K., by='region')

		tot.richness = length(unique(all.populations$spp.name[all.populations$extant == 1]))

		#if (is.element(curr.time,print.times)==T) {print(c(curr.time,nrow(all.populations),date(),tot.richness))} else{}

		## create separate subdirectory for each sim for writing output
		sim_out_dir = paste("./raw_sim_output/sim", sim, "_out", sep = "")
		if (!file.exists(sim_out_dir)) { dir.create(sim_out_dir) }
		    
    # weed out extinct populations from all.populations, but keep track of the number of extinct pops
    # and the timing of their extinction
    if (length(all.populations$extant[all.populations$extant == 0]) > 500) {

			extinct.pops = subset(all.populations, extant == 0)
			tot.extinct.pops = tot.extinct.pops + nrow(extinct.pops)
			write.csv(extinct.pops, paste(sim_out_dir, "/temp.extinct.sim.", sim, ".time.", curr.time, ".csv", sep=""), row.names=F, quote=F)
			all.populations = subset(all.populations, extant == 1)
			extinct.pops.output.times = c(extinct.pops.output.times, curr.time)

		} 

    # Stop simulation if total richness has exceeded the specified threshold
		if (tot.richness >= max.richness) {break} else{}

	} 
	# end main time loop ---------------------------------------------------------------------

  
	# Preparing simulation output
  # ----------------------------------------------------------------------------
  colnames(time.richness) = c("time", "region", "spp.rich", "reg.env")
	time.richness = as.data.frame(time.richness)

	## collect all extinct populations and combine with the current all.populations
	correct.num.rows = tot.extinct.pops + nrow(all.populations)

	if (length(extinct.pops.output.times) > 0) {
		all.pops.row.id = 1
		all.pops.out = as.data.frame(matrix(-999, ncol = ncol(all.populations), nrow = correct.num.rows))
		colnames(all.pops.out) = colnames(all.populations)
		all.pops.out[all.pops.row.id:nrow(all.populations), ] = all.populations
		all.pops.row.id = all.pops.row.id + nrow(all.populations)

		for (out.times in extinct.pops.output.times) {
		  extinct.in = read.csv(paste(sim_out_dir, "/temp.extinct.sim.",sim,".time.", as.integer(out.times),".csv",sep=""),header=T)
			all.pops.out[all.pops.row.id:(all.pops.row.id + nrow(extinct.in) - 1), ] = extinct.in
			all.pops.row.id = all.pops.row.id + nrow(extinct.in)
			rm('extinct.in')
		}
		# Remove temporary files with extinct populations
		files = list.files(sim_out_dir)
		temp.files = files[grep(paste("temp.extinct.sim.", sim, sep = ""), files)]
		file.remove(paste(sim_out_dir, "/", temp.files, sep = ""))
		
    
		all.populations = all.pops.out
		spp.extinct.times = as.matrix(tapply(all.populations$time.of.sp.extinction, all.populations$spp.name, FUN = min)) 
		spp.extinct.times = cbind(as.numeric(rownames(spp.extinct.times)), spp.extinct.times[, 1]) 
    colnames(spp.extinct.times) = c("spp.name", "time.of.sp.extinction")
		all.populations = all.populations[ , -which(colnames(all.populations) == "time.of.sp.extinction")]
		all.populations = merge(all.populations, spp.extinct.times, by = "spp.name")

	} # end if 

	#print(c('this should be zero', correct.num.rows - nrow(subset(all.populations, spp.name != -999)))) # this should be zero

	# exclude the two extreme spatial bins which may suffer boundary effects (regions 0 and 11 for our implementation)
  all.populations = subset(all.populations, region %in% 1:(num.of.bins - 1))

	## write outputs
	phylo.out = make.phylo.jimmy.fun(t = curr.time, edge.length.out = edge.length, edge.out = edge ,stem.depth.out = stem.depth )
	write.csv(all.populations, paste(sim_out_dir, "/SENC_all.pops_sim", sim, ".csv", sep=""), quote = F, row.names = F)
	write.csv(time.richness, paste(sim_out_dir, "/SENC_time.rich_sim", sim, ".csv", sep=""), quote = F, row.names = F)
	write.tree(phylo.out, paste(sim_out_dir, "/SENC_phylo_sim", sim, ".tre", sep = ""))
	end.params = data.frame(sim.id = sim,
                          status = 'completed',
                          reg.of.origin = region.of.origin,
                          w= w,
                          alpha = alpha,
                          beta = beta,
                          sigma_E = sigma_E,
                          carry.cap = carry.cap,
                          energy.gradient = energy.gradient,
                          max.K = max.K,
                          num.of.bins = num.of.bins,
                          max.time = max.time,
                          max.richness = max.richness,
                          temperate_disturb_intensity = temperate_disturb_intensity,
                          tropical_disturb_intensity = tropical_disturb_intensity,
                          disturb_frequency = disturb_frequency,
                          specn.gradient = specn.gradient,
                          specn.factor = specn.factor) 
	write.csv(end.params, paste(sim_out_dir, "/SENC_params.out_sim", sim, ".csv", sep = ""), quote = F, row.names = F)

	return(list(all.populations = all.populations, 
              time.richness = time.richness,
              phylo.out = phylo.out,
              params.out = end.params))

} # end senc_sim_fun

