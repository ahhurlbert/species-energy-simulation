## Simulation function

senc_sim_fun = function(sim.matrix, sim) {

  # Assign simulation parameters from the sim.matrix based on the sim.id
  # --------------------------------------------------------------------
  region.of.origin = sim.matrix$reg.of.origin[sim.matrix$sim.id == sim]  #'temperate' or 'tropical'
  
  w = sim.matrix$w[sim.matrix$sim.id == sim]  #strength of environmental filtering
  
  alpha = sim.matrix$alpha[sim.matrix$sim.id == sim]  #per capita speciation rate
  
  beta = sim.matrix$beta[sim.matrix$sim.id == sim]  #per capita dispersal rate to adjacent region
  
  sigma_E = sim.matrix$sigma_E[sim.matrix$sim.id == sim]  #strength of niche conservatism 
  
  carry.cap = sim.matrix$carry.cap[sim.matrix$sim.id == sim]  #'on' if limits exist to community abundance; else 'off'
  
  energy.gradient = sim.matrix$energy.gradient[sim.matrix$sim.id == sim]  #'on' if maximum abundance varies across the gradient; else 'off'
  
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
  
  	## additional parameters

	# abiotic environments and carrying capacities
	min.env = 0
	max.env = 40
	
  # regional environments and carrying capacities with and without an energy gradient
  reg.E.K. = data.frame(region = 0:num.of.bins, reg.env = seq(min.env, max.env, length = num.of.bins + 1))
	if (energy.gradient == 'on') { 
    reg.E.K.$carry.cap = seq(min.K, max.K, length=num.of.bins+1)
  } 
	if (energy.gradient == 'off') { 
    reg.E.K.$carry.cap = rep(max.K, num.of.bins+1)
  } 

	# environmental optimum for ancestor species and initial trait value in the ancestral environment
	all.populations = rbind(all.populations,c(reg.of.origin,1,1,reg.E.K.$reg.env[reg.E.K.$region == reg.of.origin],0,max.time+1,0,max.time+1)) 
	all.populations = merge(all.populations,reg.E.K.,by='region')
	all.populations = as.data.frame(all.populations)

	min.uni = -25 max.uni=25 # the min and max values for the random uniform distribution when there is no niche conservatism

	gamma = 0.1 # rate of decrease in extinction probability as population size increases

	total.mutations = 0 # the cumulative number of mutations

	## end additional parameters

	## start toggle-dependent functions

	# mutation functions, which depend on the presence or absence of niche conservatism
	mut.opt.fun = function(sigma_E.param, min.uni.param, max.uni.param, num.to.generate) { 
    rnorm(num.to.generate, mean = 0, sd = sigma_E.param) 
	}
	
	## end toggle-dependent functions

	## start individuals-dependent extinction, dispersal, and mutation probability functions

	extinction.fun = function (gamma.param,population.size) { exp(-gamma.param*population.size) }
	dispersal.fun = function (population.size,beta.param) { apply(as.matrix(population.size), 1, function(p) sum(rbinom(p, 1, beta.param))) }
	mutation.fun = function (population.size,alpha.param) { apply(as.matrix(population.size), 1, function(p) sum(rbinom(p, 1, alpha.param))) }

	## end individuals-dependent extinction, dispersal, and mutation probability functions

	## environmental fit function used in determining population sizes
	env.fit.fun = function (w.param,E,spp.opt) { exp(-(spp.opt - E)^2/(w.param^2)) }

	#### start main time loop

	for (curr.time in 1:max.time) {

		## update population sizes
		all.populations$env.fit = env.fit.fun(w.param=w,E=all.populations$reg.env,spp.opt=all.populations$env.opt)
		all.populations$pop.size = all.populations$env.fit*all.populations$carry.cap*all.populations$extant
		reg.pop.sizes = as.data.frame(tapply(all.populations$pop.size,all.populations$region,FUN=sum))
		reg.pop.sizes = cbind(rownames(reg.pop.sizes),reg.pop.sizes) colnames(reg.pop.sizes) = c("region","reg.individuals")
		all.populations = merge(all.populations,reg.pop.sizes,by='region')
		all.populations$pop.multiplier = all.populations$carry.cap/all.populations$reg.individuals
		all.populations$pop.multiplier[all.populations$pop.multiplier > 1] = 1
		if (carry.cap == 'off') {all.populations$pop.multiplier = 1} else{} # if there is no carrying capacity, then population size is not influenced by the actual community-level number of individuals
		all.populations$pop.size = round( all.populations$pop.size * all.populations$pop.multiplier ,digits=0)
		reg.pop.sizes = as.data.frame(tapply(all.populations$pop.size,all.populations$region,FUN=sum))
		reg.pop.sizes = cbind(rownames(reg.pop.sizes),reg.pop.sizes) colnames(reg.pop.sizes) = c("region","reg.individuals.real")
		all.populations = merge(all.populations,reg.pop.sizes,by='region')
    
    if (is.element(curr.time,disturb_times)==T) {
      
      for (reg.sampled in unique(all.populations$region[all.populations$extant==1])) {
        
        reg.disturb = reg.disturb.intensity$intensity[reg.disturb.intensity$region==reg.sampled]
        all.populations$pop.size[all.populations$region==reg.sampled & all.populations$extant==1] = round((1-reg.disturb)*all.populations$pop.size[all.populations$region==reg.sampled & all.populations$extant==1],digits=0)
        
      }  
          
    }

		## evaluate mutation across all extant species
		num.of.mutations = mutation.fun(population.size=all.populations$pop.size,alpha.param=alpha)
		all.populations$mutate = num.of.mutations

		## find new trait values for mutants
		all.populations$mutant.env.opt = all.populations$env.opt + mut.opt.fun(sigma_E,min.uni,max.uni,num.to.generate=nrow(all.populations))

		## evaluate dispersal across all extant species
		num.of.dispersals = dispersal.fun(population.size=all.populations$pop.size,beta.param=beta)
		all.populations$disperse = num.of.dispersals

		## evaluate extinction across all extant species
		extinction.probs = extinction.fun(gamma.param=gamma,population.size=all.populations$pop.size)
		all.populations$extinction = apply(as.matrix(extinction.probs), 1, function(p) rbinom(1, 1, prob=p))

		# collecting richness data, including only those species that have persisted across at least one time step, so not including new arrival via dispersal or mutation
		for (i in 0:num.of.bins) {

			time.richness[time.rich.row.id,] = c(curr.time,i,length(all.populations$region[all.populations$region==i & all.populations$extinction==0]),reg.E.K.$reg.env[reg.E.K.$region==i])
			time.rich.row.id = time.rich.row.id + 1

		}

		## build attributes used in phylogeny construction
		mutated.populations = all.populations[all.populations$mutate>0,]
		# leave only one mutant from each species. Assuming there are at most two mutations per species, choose the region of mutation randomly...see the 'sample' piece of the code
		mutated.populations = mutated.populations[is.element(rownames(mutated.populations),rownames(unique(mutated.populations[,c('spp.name','env.opt')],fromLast = sample(c(FALSE,TRUE),1))))==TRUE,]

		if (nrow(mutated.populations)>0) {total.mutations = total.mutations + 1} else{}

		if (nrow(mutated.populations)>0 & total.mutations == 1) { # this is just for the first mutation event
	
			stem.depth = as.numeric(c(curr.time,curr.time))
			all.populations = rbind(all.populations,c(mutated.populations$region,next.spp,1,mutated.populations$mutant.env.opt,curr.time,max.time+1,curr.time,max.time+1,rep(-999,(ncol(all.populations)-8))))

		} else{}

		if (nrow(mutated.populations)>0 & total.mutations > 1) { # this handles all mutations after the first mutation

			for (i in 1:nrow(mutated.populations)) {

				next.spp = next.spp + 1
				mutated.sp = mutated.populations$spp.name[i] # this is the true species name, not a node name, so it doesn't change.
				mutated.node = edge$to.node[edge$spp==mutated.sp]
				e = edge[which(edge$alive==1), ]
				edge$alive[which(edge$to.node==mutated.node)] = 0 
				edge$spp[edge$spp==mutated.sp] = -999
				edge = rbind(edge, c(mutated.node, next.node,1,mutated.sp), c(mutated.node, next.node + 1,1,next.spp))
				next.node = next.node + 2

				stem.depth <- c((stem.depth), curr.time, curr.time)
				x <- which(edge[, 2] == mutated.node)
				edge.length[x] <- curr.time - stem.depth[x]
				edge.length <- c(edge.length, NA, NA)

				all.populations = rbind(all.populations,c(mutated.populations$region[i],next.spp,1,mutated.populations$mutant.env.opt[i],curr.time,max.time+1,curr.time,max.time+1,rep(-999,(ncol(all.populations)-8))))

			}

		} else{}

		## disperse species that were selected to do so	
		disp.spp = all.populations[all.populations$disperse>0,]
		if ( nrow(disp.spp)>0 ) {
			disp.spp$region = disp.spp$region + sample(c(-1,1),length(disp.spp$region),replace=T,prob=c(0.5,0.5))
			disp.spp$time.of.origin = curr.time
			disp.spp = disp.spp[which(is.element(apply(disp.spp[,1:4],1 ,paste,collapse=","),apply(all.populations[,1:4],1 ,paste,collapse=","))==F),]
			all.populations = rbind(all.populations,disp.spp)
			all.populations = all.populations[is.element(all.populations$region,seq(0,num.of.bins,1))==TRUE,] # take out any populations that dispersed beyond the bins being evaluated
		} else{}

		## remove those species that are going extinct
		all.populations$time.of.extinction[all.populations$extinction==1 & all.populations$extant==1] = curr.time 
		extinct.this.time = unique(all.populations$spp.name[all.populations$extinction==1 & all.populations$extant==1]) # a vector of species with populations going extinct this time step
		all.populations$extant[all.populations$extinction==1 & all.populations$extant==1] = 0

		## find those species that went globally extinct this time step
		extant.spp = unique(all.populations$spp.name[all.populations$extant==1])
		global.extinct.this.step = extinct.this.time[is.element(extinct.this.time,extant.spp)==F]
		all.populations$time.of.sp.extinction[which(is.element(all.populations$spp.name,global.extinct.this.step)==T)] = curr.time

		## change the 'alive' status for those species that went globally extinct this time step
		edge$alive[is.element(edge$spp,global.extinct.this.step)==T] = 0
		edge.length[which(is.element(edge$spp,global.extinct.this.step)==T)] = curr.time - stem.depth[which(is.element(edge$spp,global.extinct.this.step)==T)]

		## update columns for carry.cap, reg.env, and reg.richness
		all.populations = all.populations[,1:8]
		all.populations = merge(all.populations,reg.E.K.,by='region')

		tot.richness = length(unique(all.populations$spp.name[all.populations$extant==1]))

		#if (is.element(curr.time,print.times)==T) {print(c(curr.time,nrow(all.populations),date(),tot.richness))} else{}

		if (length(all.populations$extant[all.populations$extant==0]) > 500) {

			extinct.pops = subset(all.populations,extant==0)
			tot.extinct.pops = tot.extinct.pops + nrow(extinct.pops)
			write.csv(extinct.pops,paste("temp.extinct.sim.",sim,".time.",curr.time,".csv",sep=""),row.names=F,quote=F)
			all.populations = subset(all.populations,extant==1)
			extinct.pops.output.times = c(extinct.pops.output.times,curr.time)

		} else{}

		if (tot.richness >= max.richness) {break} else{}

	} 

	#### end main time loop

	colnames(time.richness) = c('time','region','spp.rich','reg.env')
	time.richness = as.data.frame(time.richness)

	## collect all extinct populations and combine with the current all.populations

	correct.num.rows = tot.extinct.pops + nrow(all.populations)

	if (length(extinct.pops.output.times)>0) {

		all.pops.row.id = 1
		all.pops.out = as.data.frame(matrix(c(-999),ncol=10,nrow=correct.num.rows))
		colnames(all.pops.out) = colnames(all.populations)
		all.pops.out[all.pops.row.id:nrow(all.populations),] = all.populations
		all.pops.row.id = all.pops.row.id + nrow(all.populations)

		for (out.times in extinct.pops.output.times) {

			print(c(all.pops.row.id,date()))
			extinct.in = read.csv(paste("temp.extinct.sim.",sim,".time.",out.times,".csv",sep=""),header=T)
			all.pops.out[all.pops.row.id:(all.pops.row.id + nrow(extinct.in) - 1),] = extinct.in
			all.pops.row.id = all.pops.row.id + nrow(extinct.in)
			rm('extinct.in')

		}

		all.populations = all.pops.out
		spp.extinct.times = as.matrix(tapply(all.populations$time.of.sp.extinction,all.populations$spp.name,FUN=min)) 
		spp.extinct.times = cbind(as.numeric(rownames(spp.extinct.times)),spp.extinct.times[,1]) colnames(spp.extinct.times) = c('spp.name','time.of.sp.extinction')
		all.populations = all.populations[,-which(colnames(all.populations)=='time.of.sp.extinction')]
		all.populations = merge(all.populations,spp.extinct.times,by='spp.name')

	} else{}

	print(c('this should be zero',correct.num.rows - nrow(subset(all.populations,spp.name != -999)))) # this should be zero

	all.populations = subset(all.populations,region %in% 1:(num.of.bins-1))

	## write outputs

	phylo.out = make.phylo.jimmy.fun(t = curr.time, edge.length.out = edge.length, edge.out = edge ,stem.depth.out = stem.depth )
	write.csv(all.populations,paste("SENC_all.pops_sim",sim,".csv",sep=""),quote=F,row.names=F)
	write.csv(time.richness,paste("SENC_time.rich_sim",sim,".csv",sep=""),quote=F,row.names=F)
	write.tree(phylo.out,paste("SENC_phylo_sim",sim,".tre",sep=""))
	end.params = data.frame(sim.id=sim,status='completed',reg.of.origin=region.of.origin,w=w,alpha=alpha,beta=beta,sigma_E=sigma_E,carry.cap=carry.cap,
                          energy.gradient=energy.gradient,max.K=max.K,num.of.bins=num.of.bins,max.time=max.time,max.richness=max.richness,temperate_disturb_intensity = temperate_disturb_intensity,
                          tropical_disturb_intensity = tropical_disturb_intensity,disturb_frequency=disturb_frequency) 
	write.csv(end.params,paste("SENC_params.out_sim",sim,".csv",sep=""),quote=F,row.names=F)

	return(list(all.populations=all.populations,time.richness=time.richness,phylo.out=phylo.out,params.out=end.params))


} # end senc_sim_fun

