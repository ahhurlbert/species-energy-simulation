# Function for checking simulation output. 
# Specifically, checking to see that the total number of unique species in all.populations
# is the same as the number of species in phylo.out.
# The output contains a row for each sim with the two measures of richness and their difference.

file_dir = 'C:/SENCoutput'


region.phylo.plot = function(file_dir, output_dir) {
  require(ape)
  
  files = list.files(file_dir)
  allpops.files = files[grep('all.pops',files)]
  #phylo.files = files[grep('phylo',files)]
  
  sims = as.vector( sapply(allpops.files, function(x)
                      as.numeric(strsplit(strsplit(x, "sim")[[1]][2] , "\\.csv")[[1]][1])) )
  
  for (i in sims) {
      
    pops = read.csv(paste(file_dir,'/SENC_all.pops_sim',i,'.csv',sep=''), header=T)
    phy = read.tree(paste(file_dir,'/SENC_phylo_sim',i,'.tre',sep=''))
    params = read.csv(paste(file_dir,'/SENC_params.out_sim',i,'.csv',sep=''), header=T)
    
    if (nrow(pops) > 0) {
      
      png(paste(output_dir,'/reg_phylo_fig_sim',i,'.png',sep=''), height = 6, width = 8, units="in", res=120)  
      
      # Based on the region of origin for the simulation, determine whether to 
      # keep track of the minimum or maximum region reached by a species
      if (params$reg.of.origin=="tropical") { fun1 = function(x) { min(x,na.rm=T) } }
      if (params$reg.of.origin=="temperate") { fun1 = function(x) { max(x,na.rm=T) } }
      
      # Identify the most extreme region a given species has been found in
      sp.reg = aggregate(pops$region, by=list(pops$spp.name), function(x) fun1(x))
      names(sp.reg) = c('spp.name','region')
      
      # Regional colors
      reg.cols = data.frame(region = 1:10 , cols = c('purple','darkblue','blue','skyblue',
                                                     'green','lightgreen','yellow','orange',
                                                     'red','darkred'))
      # Merge colors to species
      sp.reg2 = merge(sp.reg, reg.cols, by = 'region', all.x=T)

      # Add in species in phylogeny but not pops
      missing.spp = phy$tip.label[!phy$tip.label %in% sp.reg$spp.name]
      sp.reg3 = rbind(sp.reg2, data.frame(region = rep(NA,length(missing.spp)), 
                                          spp.name = missing.spp,
                                          cols = rep('gray70', length(missing.spp))))
      sp.reg3$cols = as.character(sp.reg3$cols)
      
      # Plot phylogeny
      plot(phy, type="fan", show.tip.label=F)
      pch = tiplabels(pch=16,col = sp.reg3$cols[sapply(phy$tip.label, function(x) which(as.numeric(x)==sp.reg3$spp.name))])

      # Annotate plot
      if (params$carry.cap=='on' & params$energy.gradient=='on') {
        K.text = 'K gradient present'
      } else if (params$carry.cap=='on' & params$energy.gradient=='off') {
        K.text = 'K constant across regions'
      } else if (params$carry.cap=='off') {
        K.text = 'no K'
      }
      mtext(paste('Sim',params[1,1],', Origin =',params[1,3],', w =',params[1,4],', sigma =',params[1,7],
                  ',\ndisp = ',params[1,6],', specn =',params[1,5],',',K.text), 3 , cex = 1.2)
      legend('topright',legend=reg.cols$region , pch=16, col=as.character(reg.cols$cols), cex = .75)
      extant.spp = unique(pops[pops$extant==1,'spp.name'])
      mtext(paste(length(phy$tip.label),"species total;",
                  length(extant.spp),"extant;",
                  length(missing.spp),"missing"),1)
      dev.off()
      } # end if - data check
    print(paste(i,date(),sep="   "))
  } # end sim loop
} # end function
