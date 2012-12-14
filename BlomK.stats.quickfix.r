# Function for calculating Blomberg's K from the simulation files all.pops and phylo.out.
# Note that Blomberg's K is calculated as part of the stats output based on 
# analysis_workflow_Olympus.r as of 12/12/12.

# This script is being used to go through files retroactively in the event that K wasn't 
# calculated before.

file_dir = 'C:/SENCoutput'
analysis_dir = "//bioark.bio.unc.edu/hurlbertallen/manuscripts/cladevscommunity/analyses"
rundate = Sys.Date()
outfileName = paste(analysis_dir,'/BlombergK.output_',rundate,'.csv',sep='')

sim.matrix = read.csv('C:/Documents and Settings/Hurlbert/species-energy-simulation/SENC_Master_Simulation_Matrix.csv', header=T)
sim.matrix2 = sim.matrix[!(sim.matrix$carry.cap=='off' & 
                           sim.matrix$energy.gradient=='on') & 
                           sim.matrix$sim.id!=0,]

which.sims = sim.matrix2$sim.id

BlomK.stats = function(which.sims, file_dir, outfileName) {
  require(ape)
  require(phytools)
  
  unlink(outfileName)
  fileConn = file(outfileName, open = "wt")
  on.exit(close(fileConn))
  
  # Write header line for output file
  cat("sim", "BK.env", "BK.reg","\n", sep=",", file=fileConn)
  
  # Check that files exist for all of the sim #s requested; only use sims with files
  files = list.files(file_dir)
  allpops.files = files[grep('all.pops',files)]
  f.sims = as.vector( sapply(allpops.files, function(x)
    as.numeric(strsplit(strsplit(x, "sim")[[1]][2] , "\\.csv")[[1]][1])) )
  sims = which.sims[which.sims %in% f.sims]
  
  for (i in sims) {

    pops = read.csv(paste(file_dir,'/SENC_all.pops_sim',i,'.csv',sep=''), header=T)
    phy = read.tree(paste(file_dir,'/SENC_phylo_sim',i,'.tre',sep=''))
    extant.spp = unique(pops$spp.name[pops$extant==1])
    if (length(extant.spp) > 2) {
    
      # sub.pops contains only extant species
      sub.pops = subset(pops, extant==1)
      
      # drop extinct species from phy (note: drop.tip is much faster than prune.extinct.taxa)
      tips.to.drop = as.character(phy$tip.label[!phy$tip.label %in% as.character(sub.pops$spp.name)]);
      phy2 = drop.tip( phy , tips.to.drop );
      
      # Change any tips of branch length 0 to branch length 1; otherwise phylosig will bonk
      phy2$edge.length[phy2$edge.length==0] = 1
      
      # Calculate Blomberg's K for two traits: environmental optimum, and mean region of occurrence
      spp.traits = aggregate(sub.pops$region, by = list(sub.pops$spp.name, sub.pops$env.opt),
                             function(x) mean(x, na.rm=T))
      names(spp.traits) = c('spp.name','env.opt','region')
        
      spp.env = spp.traits$env.opt
      names(spp.env) = spp.traits$spp.name
      BK.env = phylosig(phy2, spp.env[phy2$tip.label], method="K")
        
      spp.reg = spp.traits$region
      names(spp.reg) = spp.traits$spp.name
      BK.reg = phylosig(phy2, spp.reg[phy2$tip.label], method="K")
        
      cat(i, BK.env, BK.reg, "\n", sep=',', file = fileConn)
      print(paste(i, date(), sep="   "))
    } # end if
  } # end sim loop
}
