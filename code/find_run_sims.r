### find sims that have run through the analysis pipeline

setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/sims.out.130204");
foo = list.files(pattern='SENC_Stats_sim'); head(foo); length(foo);
foo2 = grep('time5459',foo);
150 - length(foo2); # number of sims yet to analyze

### find sims that still need data analysis
### to be used when the analysis workflow does not finish

setwd("//constance/people/steg815/senc.analysis");
foo = list.files(pattern='NEW_Stats_sim'); head(foo); length(foo);
foo3 = sub('NEW_Stats_sim','',foo); head(foo3);
sims.analyzed = as.numeric(sub('_time_seq_root_only.csv','',foo3)); head(sims.analyzed);

all.sims = 3465:3664;

to.analyze = all.sims[which(is.element(all.sims,sims.analyzed)==F)];
to.analyze 
#write.csv(to.analyze,"sims.to.analyze.csv",quote=F,row.names=F);

### confirm all files are in specific folder

setwd("//olympus/steg815/senc.output");

foo = list.files(pattern='SENC_params.out_sim'); head(foo); length(foo);

### older code below here

setwd("C:/Users/steg815/Desktop/Stegen_PNNL/Spp-Energy-Niche-Conserv/SENC_Olympus/full-run-nov2012");


foo = list.files(); head(foo);
foo2 = foo[grep('SENC_all.pops_sim',foo)]; head(foo2);
foo3 = sub('SENC_all.pops_sim','',foo2); head(foo3);
foo4 = sub('.csv','',foo3); head(foo4);

sims.run = as.numeric(foo4);
all.sims = 1:2592;
to.run = all.sims[which(is.element(all.sims,sims.run)==F)];
write.table(t(to.run),"sims.to.run.txt",sep=" ",quote=F,row.names=F,col.names=F);
length(to.run)


setwd("//olympus/steg815/senc.output");

setwd("//constance/people/steg815/senc.analysis.out");

foo = list.files(pattern='SENC_Stats_sim'); head(foo); length(foo);
foo3 = sub('SENC_Stats_sim','',foo); head(foo3);
foo4 = sub('.csv','',foo3); head(foo4);

sims.run = as.numeric(foo4);
all.sims = 3465:4064;
to.run = all.sims[which(is.element(all.sims,sims.run)==F)];
length(to.run); head(to.run);

not.run.out = read.table('sim.3782.out'); not.run.out;

#write.table(t(to.run),"//olympus/steg815/SENC/sims.to.zip.txt",sep=" ",quote=F,row.names=F,col.names=F);
#write.table(t(2:2592),"sims.to.zip.txt",sep=" ",quote=F,row.names=F,col.names=F);

