make.phylo.jimmy.fun = function(t,edge.length.out,edge.out,stem.depth.out) {

edge.length.out[which(edge.out$alive==1)] = t - stem.depth.out[which(edge.out$alive==1)]
edge.out = as.matrix(edge.out);
    n = -1
    for (i in 1:max(edge.out[,c('from.node','to.node')])) {
        if (any(edge.out[, 'from.node'] == i)) {
            edge.out[which(edge.out[, 'from.node'] == i), 'from.node'] = n
            edge.out[which(edge.out[, 'to.node'] == i), 'to.node'] = n
            n = n - 1
        }
    }

edge.out[which(edge.out[,c('from.node','to.node')] > 0)] = 1:length(edge.out[which(edge.out[,c('from.node','to.node')] > 0)]); edge.out;
tip.label = edge.out[,'spp'][edge.out[,'spp']!=-999]; tip.label;
edge.only = edge.out[,c('from.node','to.node')]
mode(edge.only) = "character"
mode(tip.label) = "character"
obj = list(edge = edge.only, edge.length = edge.length.out, tip.label = tip.label)
class(obj) = "phylo"
phylo.out = old2new.phylo(obj)
phylo.out = read.tree(text = write.tree(phylo.out))
return(phylo.out)

};