##############################################
# Author: Allen Hurlbert
# Function: Script for summarizing simulation output, primarily to see whether certain 
#   parameter combinations were resulting in more skipped clades, fewer output rows,
#   species, or number of regions in the analysis
#

user = 'james';

if (user == 'james') {output = read.csv('//olympus/steg815/senc.output/sim.matrix.output_2012-12-14.csv',header=T); head(output);};
if (user == 'allen') {output = read.csv('//Bioark.bio.unc.edu/hurlbertallen/Manuscripts/CladeVsCommunity/analyses/sim.matrix.output_2012-12-14.csv',header=T)};

# Get rid of initial NA row
output = output[-1,]
output = droplevels(output)

# drop unwanted parameter values
output = subset(output,output$w > 1); head(output);
output = subset(output,output$alpha == 1.00E-06); head(output);
output = subset(output,output$sigma_E == 5); head(output);

output$beta_over_alpha = output$beta/output$alpha
output$skipped.clades.pct = output$skipped.clades/output$extant.S

boxwidth = .3

if (user == 'james') {pdf(paste('//olympus/steg815/senc.output/simulation_output_summary_',gsub(":","_",gsub(" ","-",Sys.time())),'.pdf',sep=''),height = 6, width = 8)};
if (user == 'allen') {pdf(paste('//Bioark.bio.unc.edu/hurlbertallen/Manuscripts/CladeVsCommunity/analyses/simulation_output_summary_',Sys.Date(),'.pdf',sep=''),height = 6, width = 8)};

par(mfrow=c(2,3),oma = c(4,1,1,1),mar = c(3,5,1,1),cex.lab = 1.5)
# Carrying capacity on or off
boxplot(output$n.regions ~ output$carry.cap, ylab = "Number of regions",boxwex= boxwidth, notch=T,col='darkblue')

boxplot(log10(output$extant.S) ~ output$carry.cap, ylab = "Log10 Number of extant species", notch=T,boxwex = boxwidth,col = 'orange2')

boxplot(log10(output$extinct.S) ~ output$carry.cap,ylab = "Log10 Number of extinct species", notch=T,boxwex = boxwidth, col = 'yellow4')

boxplot(output$skipped.clades ~ output$carry.cap, notch=T, boxwex = boxwidth, ylab = "Number of skipped clades", col='skyblue')

boxplot(output$output.rows ~ output$carry.cap, notch=T, boxwex = boxwidth, ylab = "Number of output rows", col = 'darkred')

plot(1,1,type="n",axes=F,ylab="",xlab="")
mtext("Carrying Capacity",1,outer=T, cex=1.5)

xs = c(4:7,20)
ys = c(14,15,16,21,19)

for (i in xs) {
  for (y in ys) {
    if (i %in% c(5,6,20)) { xvar = log10(output[,i])}
    else { xvar = output[,i] }
    plot(output[output$carry.cap=="off",y] ~ jitter(xvar[output$carry.cap=="off"]), ylab = names(output)[y], 
         xlab="", col='darkblue',pch=2, cex=1.5, ylim = range(output[,y],na.rm=T))
    points(output[output$carry.cap=="on",y] ~ jitter(xvar[output$carry.cap=="on"]), 
           xlab="",ylab = names(output)[y], col='red',pch=16)
  }
  plot(1,1,type="n",axes=F,ylab="",xlab="")
  legend("center",c("K on","K off"), pch = c(16,2),col = c('red','darkblue'), cex=2)
  mtext(names(output)[i], 1, outer=T, cex=2)
}
dev.off()
