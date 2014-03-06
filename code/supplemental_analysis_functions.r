# Author: Allen Hurlbert

# Functions used in the analysis of simulation output based on Hurlbert & Stegen 2014, Ecology Letters


#Function for unzipping output files, and storing them in a list
# --output_dir is the directory where the zipped output files reside

output.unzip = function(output_dir, sim_ID) {
  require(ape)
  
  # check that raw simulation output files exist
  # file.access() returns value of 0 when mode=0 if file exists
  if (file.access(paste(output_dir, '/sim', sim_ID, '_out/SENC_all.pops_sim', sim_ID, '.csv', sep = ''), mode = 0) == 0) { 
    
    all.populations = read.csv(paste(output_dir, '/sim', sim_ID, '_out/SENC_all.pops_sim', sim_ID, '.csv', sep = ''), header=T)
    params.out = read.csv(paste(output_dir, '/sim', sim_ID, '_out/SENC_params.out_sim', sim_ID, '.csv', sep = ''), header=T)
    phylo.out = read.tree(paste(output_dir, '/sim', sim_ID, '_out/SENC_phylo_sim', sim_ID, '.tre', sep = ''))
    time.richness = read.csv(paste(output_dir, '/sim', sim_ID, '_out/SENC_time.rich_sim', sim_ID, '.csv', sep = ''), header=T)
    sim.results = list(all.populations=all.populations, 
                       params.out=params.out, 
                       phylo.out=phylo.out, 
                       time.richness=time.richness)
  } else if( file.access(paste(output_dir,'/senc.out.',sim_ID,'.zip',sep=''), mode=0) == 0 ) {
    #if not, then what about zipped sim output?
    unzipped.filenames = unzip(paste(output_dir,'/senc.out.',sim_ID,'.zip',sep=''), 
                               exdir = paste(output_dir, '/sim', sim_ID, '_out', sep = ''))
    all.populations = read.csv(unzipped.filenames[1], header=T)
    params.out = read.csv(unzipped.filenames[2],header=T)
    phylo.out = read.tree(unzipped.filenames[3])
    time.richness = read.csv(unzipped.filenames[4], header=T)
    sim.results = list(all.populations=all.populations, 
                       params.out=params.out, 
                       phylo.out=phylo.out, 
                       time.richness=time.richness)
  } else {
    print("No simulation output exists for that sim ID in the directory specified")
    sim.results = NULL
  }
    return(sim.results)
}

############################################################################################################

# Function for making a phylogeny (of class 'phylo') from the simulation output

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

###########################################################################################################
  
# For all possible subclades
  # Over all available time slices
    # In each region                                    [function = regional.calc()]
      # Function inputs: (1) a dataframe with site (region), species, and time of origin columns, 
      #                  (2) a 'phylo' class object for those species, and (3) the amount of time
      #                  simulations were run for ('max.time'), and
      # Calculate richness
      # Calculate MRD
        # This part of the function calculates the mean root distance for each assemblage in each site.
        # Modified from code originally by ELIOT MILLER 25 AUGUST 2011
      # Calculate PSV
        # Helmus et al. 2007's Phylogenetic Species Variability metric
        # based on the phylogenetic variance/covariance matrix
        # While Helmus et al.'s equation is as follows
        # PSV = (n*trace(C) - sum(C))/(n*(n-1)) 
        # this only holds for a tree with contemporaneous tips where C is the correlation matrix.
        # When the tree is not ultra-metric, Algar et al. 2009 (Ecol Lett) show the calculation
        # using the variance-covariance matrix, V:
        # PSV = (n*trace(V) - sum(V))/(trace(V)*(n-1))
      # Calculate ...?
    # Across regions calculate                          [function = xregion.analysis()]
      # Calculate environment-richness slope
      # Calculate time-richness slope
      # Calculate environment-MRD slope
      # Calculate environment-PSV slope
      # Calcualte MRD-richness slope
      # Calculate PSV-richness slope
      # Calculate partial correlation btw environment-richness after accounting for MRD and PSV (?)
regional.calc = function(sub.populations, phylo.out, max.time)
{
  #Required libraries
  #require(plyr)
  #require(picante)
  #require(ape)
  
  if(class(phylo.out) != "phylo") {
    stop("Second argument must be of class 'phylo'\n")
  }
  if(class(sub.populations) != "data.frame") {
    stop("First argument must be a dataframe\n")
  } 
  if(class(sub.populations) == "data.frame" & ncol(sub.populations) != 5) {
    stop("First argumument should have 5 cols: sites, species, time of origin, and reg.env\n")
  }
  if(class(max.time) != "integer") {
    stop("max.time must be an integer\n")
  }
  
  global.clade.richness = length(unique(sub.populations$spp.name))
  clade.extant.richness = length(unique(sub.populations[sub.populations$extant==1,'spp.name']))
  
  # Calculate the time of origin of the focal clade within each region
  overall.origin.time = phylo.out$origin.time;
  origin.by.region = aggregate(sub.populations$time.of.origin, by=list(sub.populations$region), min)
  names(origin.by.region) = c('region','clade.origin.time')
  origin.by.region$clade.origin.time[origin.by.region$clade.origin.time < overall.origin.time] = overall.origin.time
  reg.time = data.frame(region = origin.by.region$region, 
                        time.in.region = max.time - origin.by.region$clade.origin.time)
  
  #I added this line to remove duplicate instances of same species due to repeated colonizations
  sub.pops = unique(sub.populations[,c('region','spp.name','reg.env')])
  
  #MRD
  phylo.bl1 <- compute.brlen(phylo.out, 1)
  all.dist <- dist.nodes(phylo.bl1)
  root.dist <- all.dist[length(phylo.out$tip.label)+1, 1:length(phylo.out$tip.label)]
  tips.to.root <- data.frame(spp.name=phylo.out$tip.label,root.dist)
  MRD.ini <- merge(sub.pops, tips.to.root, sort = FALSE)
  MRD.ini <- MRD.ini[order(MRD.ini$region), ]
  MRD.table <- ddply(MRD.ini, "region", summarise, RD=sum(root.dist), richness=length(unique(spp.name)))
  MRD <- data.frame(region=MRD.table$region, richness=MRD.table$richness, MRD=MRD.table$RD/MRD.table$richness)
  
  MRD2 <- merge(MRD, reg.time, by = "region", all = T)
  
  #PSV
  Vmatrix = vcv(phylo.out, corr=F)
  psvs = matrix(NA, ncol=2)
  #Can only calculate where number of species in a region > 1
  reg.S = data.frame(table(sub.pops$region))
  regions = as.numeric(as.character(reg.S[reg.S$Freq > 1, 'Var1']))
  for (i in regions) {
    subset.sp = sub.populations[sub.populations$region==i, 'spp.name']
    index = as.numeric(row.names(Vmatrix)) %in% subset.sp
    v.matrix = Vmatrix[index,index]
    n = nrow(v.matrix)
    psv = (n*sum(diag(v.matrix)) - sum(v.matrix))/(sum(diag(v.matrix))*(n-1))
    psvs = rbind(psvs, c(i,psv))
  }
  PSVs = data.frame(matrix(psvs[-1,],ncol=2)) #remove row of NAs
  names(PSVs) = c('region','PSV')
    
  MRD.PSV.out = merge(MRD2, PSVs, by = 'region',all=T)
  MRD.PSV.out$clade.origin = rep(overall.origin.time, nrow(MRD.PSV.out))
  MRD.PSV.out = unique(merge(MRD.PSV.out,sub.populations[,c('region','reg.env')],by='region'))
    
  MRD.PSV.out$global.richness = global.clade.richness
  MRD.PSV.out$extant.richness = clade.extant.richness
  MRD.PSV.out$max.RD = max(root.dist)
  return(MRD.PSV.out)
}
  


xregion.analysis = function(region.summary) {
  
  n.regions = nrow(region.summary)

  corr.output = data.frame(r.time.rich = NA,  #corr betw richness and time-in-region
                           p.time.rich = NA,
                           r.env.rich = NA,   #corr betw richness and environment
                           p.env.rich = NA,
                           r.MRD.rich = NA,   #corr betw richness and MRD
                           p.MRD.rich = NA,
                           r.PSV.rich = NA,   #corr betw richness and PSV
                           p.PSV.rich = NA,
                           r.env.MRD = NA,    #corr betw environment and MRD
                           p.env.MRD = NA,   
                           r.env.PSV = NA,    #corr betw environment and PSV
                           p.env.PSV = NA,
                           r.ext.reg = NA,    #corr betw extinction rate and region
                           p.ext.reg = NA,
                           r.rich.ext = NA,   #corr betw richness and extinction rate
                           p.rich.ext = NA,
                           n.regions = n.regions,
                           clade.origin.time = region.summary$clade.origin[1],
                           global.richness = region.summary$global.richness[1],
                           extant.richness = region.summary$extant.richness[1],
                           max.RD = region.summary$max.RD[1])
    
  if(nrow(unique(region.summary[!is.na(region.summary$time.in.region) & !is.na(region.summary$richness),c('time.in.region','richness')])) > 2) {
    r1 = cor.test(region.summary$time.in.region, region.summary$richness, method = "pearson")
    corr.output$r.time.rich = r1$estimate
    corr.output$p.time.rich = r1$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$reg.env) & !is.na(region.summary$richness),c('reg.env','richness')])) > 2) {
    r2 = cor.test(region.summary$reg.env, region.summary$richness, method = "pearson")
    corr.output$r.env.rich = r2$estimate
    corr.output$p.env.rich = r2$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$MRD) & !is.na(region.summary$richness),c('MRD','richness')])) > 2) {
    r3 = cor.test(region.summary$MRD, region.summary$richness, method = "pearson")
    corr.output$r.MRD.rich = r3$estimate
    corr.output$p.MRD.rich = r3$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$PSV) & !is.na(region.summary$richness),c('PSV','richness')])) > 2) {
    r4 = cor.test(region.summary$PSV, region.summary$richness, method = "pearson")
    corr.output$r.PSV.rich = r4$estimate
    corr.output$p.PSV.rich = r4$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$reg.env) & !is.na(region.summary$MRD),c('reg.env','MRD')])) > 2) {
    r5 = cor.test(region.summary$reg.env, region.summary$MRD, method = "pearson")
    corr.output$r.env.MRD = r5$estimate
    corr.output$p.env.MRD = r5$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$reg.env) & !is.na(region.summary$PSV),c('reg.env','PSV')])) > 2) {
    r6 = cor.test(region.summary$reg.env, region.summary$PSV, method = "pearson")
    corr.output$r.env.PSV = r6$estimate
    corr.output$p.env.PSV = r6$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$extinction.rate) & !is.na(region.summary$region),c('reg.env','PSV')])) > 2) {
    r7 = cor.test(region.summary$extinction.rate, region.summary$region, method = "pearson")
    corr.output$r.ext.reg = r7$estimate
    corr.output$p.ext.reg = r7$p.value
  }
  if(nrow(unique(region.summary[!is.na(region.summary$richness) & !is.na(region.summary$extinction.rate),c('reg.env','PSV')])) > 2) {
    r8 = cor.test(region.summary$richness, region.summary$extinction.rate, method = "pearson")
    corr.output$r.rich.ext = r8$estimate
    corr.output$p.rich.ext = r8$p.value
  }
  return(corr.output)
}
  

# Testing of above functions to make sure they give the same output as that in Figure 1
# of Algar et al. 2009 (Ecology Letters 12: 57-65)

if(0){
# Tree 1 (strong niche conservatism)
akctree1 = read.tree(text="((1,2),((((4,5),3),6),((7,8),((10,11),9))));")
akctree1$edge.length = rep(1, nrow(akctree1$edge))
com1 = data.frame(rbind(cbind(rep(1,6),1:6),cbind(rep(2,3),7:9),cbind(rep(3,2),10:11)))
names(com1) = c('region','tipnames')
com1$time.of.origin = rep(1,nrow(com1))

# Tree 2 (weak niche conservatism)
akctree2 = read.tree(text="(1,((2,(3,4)),(5,6)));")
akctree2$edge.length = rep(1,nrow(akctree2$edge))
com2 = data.frame(rbind(cbind(rep(1,2),1:2),cbind(rep(2,2),3:4),cbind(rep(3,2),5:6)))
names(com2) = c('sites','spp')
com2$time.of.origin = rep(1,nrow(com2))

regional.calc(com1, akctree1, 1)
#expected output: 
# region richness MRD time.in.region       PSV clade.origin
# 1      1        6 3.5             0 0.6761905            1
# 2      2        3 4.0             0 0.4166667            1
# 3      3        2 5.0             0 0.2000000            1

regional.calc(com2, akctree2, 1)
#expected output:
# region richness MRD time.in.region       PSV clade.origin
# 1      1        2   2             0 1.0000000            1
# 2      2        2   4             0 0.2500000            1
# 3      3        2   3             0 0.3333333            1
}


###########################################################################################################################

# Function for calculating extinction rate
# This is in units of extinctions per population per time, calculated as the 
# total number of extinctions in a region divided by the total number of populations 
# that ever appeared (whether through speciation or dispersal) divided by the time
# since the first colonization of that region.
# The argument should be an all.populations dataframe subsetted down to the appropriate
# timeslice (i.e., without species with time.of.origin > t)

extinct.calc = function(all.populations, timeslice) {
  sub.populations = subset(all.populations, time.of.origin <= timeslice)
  num.all.lineages = data.frame(table(sub.populations$region))
  extinct.lineages = aggregate(sub.populations$time.of.extinction,by=list(sub.populations$region), function(x) sum(x<timeslice))
  time.col = aggregate(sub.populations$time.of.origin,by=list(sub.populations$region), min)
  output = data.frame(region = extinct.lineages$Group.1, extinct.pops=extinct.lineages$x, 
                      total.pops=num.all.lineages$Freq, time.in.region.pops=timeslice-time.col$x)
  output$extinction.rate = output$extinct.pops/output$total.pops/output$time.in.region.pops
  return(output)   
}

##########################################################################################################################

# This function is a modification of the maxlik.betasplit() function
# in the apTreeshape package for calculating beta, a metric of phylogenetic tree imbalance.

maxlik.betasplit.AH = function (phylo, up = 10, remove.outgroup = FALSE, confidence.interval = "none", 
          conf.level = 0.95, size.bootstrap = 100) 
{
  vrais.aldous.fin <- function(i, n, b) {
    # Code commented out below is from the original maxlik.betasplit function.
    # Due to underflow errors, beta was switched out for lbeta in the 
    # manner of Purvis et al. 2011, Phil Trans Roy Soc B 366: 2462-2477
    # (Purvis, pers. comm.)
    #aux <- beta(b + i + 1, b + n - i + 1)/beta(i + 1, n - 
    #                                             i + 1)
    #if (is.na(aux) | (aux == Inf)) 
    #  aux <- (i/n)^b * (1 - i/n)^(b)
    aux <- exp(lbeta(b + i + 1, b + n - i + 1) - lbeta(i + 1, n - i + 1))
    return(aux)
  }
  bbalance <- function(phylo) {
    return(t(apply(balance(phylo), FUN = function(x) {
      c(x[1], x[1] + x[2])
    }, MARGIN = 1)))
  }
  renorm.aldous <- function(n, beta) {
    return(sum(sapply(1:(n - 1), FUN = vrais.aldous.fin, 
                      n = n, b = beta)))
  }
  vrais.aldous.renorm <- function(i, n, beta) {
    return(vrais.aldous.fin(i, n, beta)/renorm.aldous(n, 
                                                      beta))
  }
  logvrais.aldous.phylo <- function(b, phylo, remove.outgroup = TRUE) {
    if (class(phylo) == "treeshape") 
      bal <- bbalance(as.phylo(phylo))
    if (class(phylo) == "phylo") 
      bal <- bbalance(phylo)
    if (remove.outgroup) {
      if ((bal[1, 1] <= 2) || ((bal[1, 2] - bal[1, 1]) <= 
                                 2)) 
        bal <- bal[-1, ]
    }
    return(sum(log(apply(bal, FUN = function(x, b) {
      return(vrais.aldous.renorm(x[1], x[2], b))
    }, b = b, MARGIN = 1))))
  }
  logvrais.aldous.bal <- function(b, bal) {
    return(sum(log(apply(bal, FUN = function(x, b) {
      return(vrais.aldous.renorm(x[1], x[2], b))
    }, b = b, MARGIN = 1))))
  }
  if (class(phylo) == "treeshape") 
    bal <- bbalance(as.phylo(phylo))
  if (class(phylo) == "phylo") 
    bal <- bbalance(phylo)
  if (class(phylo) != "phylo" && class(phylo) != "treeshape") {
    print("The phylogeny shall be of class phylo or treeshape")
    return
  }
  if (remove.outgroup) {
    if ((bal[1, 1] <= 2) || ((bal[1, 2] - bal[1, 1]) <= 2)) 
      bal <- bal[-1, ]
  }
  nb.tip <- max(bal[1, ])
  optim_lik_aldous <- function(phylo, remove.outgroup) {
    optimize(f = function(x) {
      logvrais.aldous.phylo(x, phylo, remove.outgroup)
    }, lower = -2, upper = up, maximum = TRUE)
  }
  optim_lik_aldous_bal <- function(bal) {
    optimize(f = function(x) {
      logvrais.aldous.bal(x, bal)
    }, lower = -2, upper = up, maximum = TRUE)
  }
  res <- optim_lik_aldous(phylo, remove.outgroup)
  if (confidence.interval == "bootstrap") {
    nb_b <- dim(bal)[1]
    fun_aux <- function() {
      optim_lik_aldous_bal(bal[sample(1:nb_b, size = nb_b, 
                                      replace = TRUE), ])$maximum
    }
    thebeta <- replicate(size.bootstrap, fun_aux())
    up.conf <- 1 - ((1 - conf.level)/2)
    low.conf <- (1 - conf.level)/2
    conf_interval <- quantile(thebeta, c(low.conf, up.conf))
  }
  if (confidence.interval == "profile") {
    if ((res$objective - 1.92) - logvrais.aldous.phylo(-2, 
                                                       phylo, remove.outgroup) < 0) 
      low.conf <- (-2)
    else low.conf <- uniroot(f = function(x) {
      (res$objective - 1.92) - logvrais.aldous.phylo(x, 
                                                     phylo, remove.outgroup)
    }, lower = -2, upper = res$maximum)$root
    if ((res$objective - 1.92) - logvrais.aldous.phylo(up, 
                                                       phylo, remove.outgroup) < 0) 
      up.conf <- up
    else up.conf <- uniroot(f = function(x) {
      (res$objective - 1.92) - logvrais.aldous.phylo(x, 
                                                     phylo, remove.outgroup)
    }, lower = res$maximum, upper = up)$root
    conf_interval <- c(low.conf, up.conf)
  }
  if (confidence.interval == "none") {
    conf_interval <- NULL
  }
  return(list(max_lik = res$maximum, conf_interval = conf_interval))
}