# Function for checking simulation output. 
# Specifically, checking to see that the total number of unique species in all.populations
# is the same as the number of species in phylo.out.
# The output contains a row for each sim with the two measures of richness and their difference.

file_dir = 'C:/SENCoutput'


region.envopt.phylo.plot = function(file_dir, output_dir, sims=F, extant.only=F) {
  require(ape)
  
  if (!sims) {
    files = list.files(file_dir)
    allpops.files = files[grep('all.pops',files)]
    #phylo.files = files[grep('phylo',files)]
  
    sims = as.vector( sapply(allpops.files, function(x)
                      as.numeric(strsplit(strsplit(x, "sim")[[1]][2] , "\\.csv")[[1]][1])) )
  } 
  
  for (i in sims) {
      
    pops = read.csv(paste(file_dir,'/SENC_all.pops_sim',i,'.csv',sep=''), header=T)
    phy = read.tree(paste(file_dir,'/SENC_phylo_sim',i,'.tre',sep=''))
    params = read.csv(paste(file_dir,'/SENC_params.out_sim',i,'.csv',sep=''), header=T)
    
    if(extant.only) {
      extant.ornot = aggregate(pops$extant,by=list(pops$spp.name),sum)
      extinct.species = as.character(extant.ornot[extant.ornot$x==0,'Group.1'])
      sub.species = as.character(unique(pops$spp.name));
      sub.species2 = sub.species[!sub.species %in% extinct.species]
      tips.to.drop = as.character(phy$tip.label[!phy$tip.label %in% sub.species2]);
      sub.phy = drop.tip(phy,tips.to.drop);
    }
    
    if (nrow(pops) > 0) {
      
      png(paste(output_dir,'/reg_envopt_phylo_fig_sim',i,'.png',sep=''), height = 6, width = 12, units="in", res=72)  
      par(mfrow=c(1,2), oma = c(2,0,4,0), mar = c(1,1,1,2))
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
      # Merge colors to species for regions
      sp.reg2 = merge(sp.reg, reg.cols, by = 'region', all.x=T)
      sp.reg2$cols = as.character(sp.reg2$cols)

      # Drop species not in pops from phylogeny (these are boundary species only present in regs 0,11)
      missing.spp = sub.phy$tip.label[!sub.phy$tip.label %in% sp.reg$spp.name]
      phy2 = drop.tip(sub.phy, missing.spp)
      
      # Plot phylogeny coded by region
      plot(phy2, type="fan", show.tip.label=F)
      tiplabels(pch=16,col = sp.reg2$cols[sapply(phy2$tip.label, function(x) which(as.numeric(x)==sp.reg2$spp.name))])
      legend('topright',legend = 1:10, pch = 16, col = as.character(reg.cols$cols))
      legend('topleft','Region', bty='n', cex = 1.2)
      
      # Environmental optimum colors
      envopt.cols = colorRampPalette(c('darkblue','skyblue','gray95','pink','red','darkred'))(40)
      
      spp.envopts = unique(pops[,c('spp.name','env.opt')])
      
      # Rescale environmental optimum values to range between 0 and 40 for converting to colors
      envopts = 39*(spp.envopts$env.opt - min(spp.envopts$env.opt))/(max(spp.envopts$env.opt) - min(spp.envopts$env.opt)) + 1
      names(envopts) = spp.envopts$spp.name
      
      # Plot phylogeny coded by environmental optimum
      plot(phy2, type="fan", show.tip.label=F)
      tiplabels(pch = 16, col = envopt.cols[round(envopts[phy2$tip.label])])
      legend('topright',legend = round(seq(min(spp.envopts$env.opt),max(spp.envopts$env.opt),length.out = 5),1), 
             pch = 16, col = envopt.cols[seq(1,40,length.out=5)])
      legend('topleft','Environmental\nOptimum', bty='n', cex = 1.2)
      
      # Annotate plot
      if (params$carry.cap=='on' & params$energy.gradient=='on') {
        K.text = 'K gradient present'
      } else if (params$carry.cap=='on' & params$energy.gradient=='off') {
        K.text = 'K constant across regions'
      } else if (params$carry.cap=='off') {
        K.text = 'no K'
      }
      mtext(paste('Sim',params[1,1],', Origin =',params[1,3],', w =',params[1,4],', sigma =',params[1,7],
                  ',\ndisp = ',params[1,6],', specn =',params[1,5],',',K.text), 3 , cex = 1.75, outer = T)
      
      extant.spp = unique(pops[pops$extant==1,'spp.name'])
      mtext(paste(length(phy$tip.label),"species total;",
                  length(extant.spp),"extant;",
                  length(missing.spp),"missing"),1, outer = T, cex = 1.5)
      dev.off()
      } # end if - data check
    print(paste(i,date(),sep="   "))
  } # end sim loop
} # end function
