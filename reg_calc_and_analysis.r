# Author: Allen Hurlbert
# Date: 23 October 2012
  
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

###################
###################
# NOTE: I HAVE TEMPORARILY COMMENTED OUT LINES INVOLVING overall.origin.time
#
# This still needs to be fixed and merged back into the output at the end


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
  if(class(sub.populations) == "data.frame" & ncol(sub.populations) != 4) {
    stop("First argumument should have 4 cols: sites, species, time of origin, and reg.env\n")
  }
  if(class(max.time) != "integer") {
    stop("max.time must be an integer\n")
  }
  

  #Calculate the time of origin of the focal clade within each region
  #overall.origin.time = sub.clade.phylo$origin.time;
  origin.by.region = aggregate(sub.populations$time.of.origin, by=list(sub.populations$region), min)
  names(origin.by.region) = c('region','clade.origin.time')
  #origin.by.region$clade.origin.time[origin.by.region$clade.origin.time < overall.origin.time] = overall.origin.time
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
  MRD.table <- ddply(idata.frame(MRD.ini), "region", summarise, RD=sum(root.dist), richness=length(unique(spp.name)))
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
  #MRD.PSV.out$clade.origin = rep(overall.origin.time, nrow(MRD.PSV.out))
  MRD.PSV.out = unique(merge(MRD.PSV.out,sub.populations[,c('region','reg.env')],by='region'))
  return(MRD.PSV.out)
}
  


xregion.analysis = function(region.summary) {
  
  n.regions = nrow(region.summary)

  corr.output = data.frame(r.time.rich = NA,
                           p.time.rich = NA,
                           r.env.rich = NA,
                           p.env.rich = NA,
                           r.MRD.rich = NA,
                           p.MRD.rich = NA,
                           r.PSV.rich = NA,
                           p.PSV.rich = NA,
                           r.env.MRD = NA,
                           p.env.MRD = NA,
                           r.env.PSV = NA,
                           p.env.PSV = NA,
                           n.regions = n.regions,
                           clade.origin.time = region.summary$clade.origin[1])
    
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
  return(corr.output)
}
  


#######################################################################################
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