# Script for compiling output across simulations for the root clade
# finding the root clade as the clade with the most species and the earliest origin time
# similar to compile.firstline.r, but using richness and age to make sure we get the root clade

Allen = 0;

if (Allen == 1) {

file_dir = '//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses'
filetext = "SENC_Stats_sim"

}

if (Allen == 0) {
  
  file_dir = 'C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204'
  filetext = "SENC_Stats_sim"
  
}


# get a list of the proper files, pull out the root clade of each and store together in a df

compile.firstlines = function(file_dir, filetext) {
  files = list.files(file_dir)
  stats.files = files[grep(filetext,files)]
  simstats.finaltime = c()
  for (i in stats.files) {
    stats = read.csv(paste(file_dir,'/',i,sep=''), header=T)
    if (nrow(stats) > 0) {
      
      root.from.richness = which.max(stats$extant.richness)
      root.from.age = which.min(stats$clade.origin.time)
      
      if (root.from.richness == root.from.age) { 
      
        simstats.finaltime = rbind(simstats.finaltime, stats[root.from.richness,])
        
      } else { print(c("Error: largest clade is not the oldest clade",paste("Sim = ",stats$sim[1],sep="")))  }
      
            
    } else {
      sim = as.numeric( strsplit( strsplit(i,"m")[[1]][2] , "\\.")[[1]][1] )
      simstats.finaltime = rbind(simstats.finaltime, c(sim, rep(NA,(ncol(simstats.finaltime)-1))))
    }
  }
  return(simstats.finaltime)
}

root.comp = compile.firstlines(file_dir=file_dir,filetext=filetext);
head(root.comp); dim(root.comp); range(root.comp$sim);

write.csv(root.comp,paste(file_dir,"/rootclade_stats_sims",min(root.comp$sim),"-",max(root.comp$sim),".csv",sep=""),quote=F,row.names=F)