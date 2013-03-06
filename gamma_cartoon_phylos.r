
setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/senc_ms_figs");

library('ape');

pdf("cartoon_gammas.pdf",height=5,width=15);

  par(mfrow=c(1,3))

  neg = read.csv("neg_gamma_phylo.csv",header=T,row.names=1); neg
  neg.phylo = as.phylo(hclust(as.dist(neg)));
  plot(neg.phylo,show.tip.label=F,edge.width=2)
  title(paste("Gamma = ",round(gammaStat(neg.phylo),digits=0),sep=""),cex.main=2);

  bal = read.csv("balanced_phylo.csv",header=T,row.names=1); bal
  bal.phylo = as.phylo(hclust(as.dist(bal)));
  plot(bal.phylo,show.tip.label=F,edge.width=2);
  title(paste("Gamma = ",round(gammaStat(bal.phylo),digits=0),sep=""),cex.main=2);

  pos = read.csv("pos_gamma_phylo.csv",header=T,row.names=1); pos
  pos.phylo = as.phylo(hclust(as.dist(pos)));
  plot(pos.phylo,show.tip.label=F,edge.width=2)
  title(paste("Gamma = ",round(gammaStat(pos.phylo),digits=0),sep=""),cex.main=2);

dev.off();
