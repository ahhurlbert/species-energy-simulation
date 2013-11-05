# find maximum times for specified simulations

setwd("//constance/people/steg815/senc.sims.130322");

sims = 3465:3664;

max.times = numeric();

for (sim in sims) {
  
  time.rich.temp = read.csv(paste("SENC_time.rich_sim",sim,".csv",sep=""));
  max.times = rbind(max.times,c(sim,max(time.rich.temp$time)));
  print(sim);
  
}

max.times = as.data.frame(max.times);
colnames(max.times) = c('sim','time');
head(max.times);

min(max.times$time);
hist(max.times$time);
