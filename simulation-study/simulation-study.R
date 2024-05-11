source("examples/Parameters/BDI/Equilibrium (N=100).R")

### REMINDER: Reusing estimates?

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
tau = 0.1
box = get_box("BDI",theta,tout_list,tau,x,obs_list,0.9)
lower_list = box$lower
upper_list = box$upper

pdf("simulation-study/BDI/N=100/different-m/plot.pdf")
plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
dev.off()

### M = 10
pdf("simulation-study/BDI/N=100/different-m/M=10.pdf")
M = 10
N = 10000
frac_sims = list(sims=rep(0,N),times=rep(0,N))
rb_sims = list(sims=rep(0,N),times=rep(0,N))

for(i in 1:N){
  start.time = Sys.time()
  frac_sims$sims[i] = frac_list("BDI",x, theta, tout_list, obs_list,M)
  end.time = Sys.time()
  frac_sims$times[i] = end.time-start.time
}

for(i in 1:N){
  start.time = Sys.time()
  rb_sims$sims[i] = RB_list("BDI",x, theta, tout_list, lower_list, upper_list, tau,obs_list,M)
  end.time = Sys.time()
  rb_sims$times[i] = end.time-start.time
}
eff_rb = formatC(var(rb_sims$sims)*sum(rb_sims$times),format='e')
eff_frac = formatC(var(frac_sims$sims)*sum(frac_sims$times),format='e')

p1 = hist(frac_sims$sims,breaks=25,prob=T,col=rgb(0,0,1,1/4), xlim=c(0,0.55),ylim=c(0,40),
main='Estimates of likelihood',cex.lab=1.5,cex.main=1.5,cex.axis=1.5,xlab='likelihood',ylab='density')
p2 = hist(rb_sims$sims,breaks=25,freq=F,prob=T,add=T, col=rgb(1,0,0,1/4), xlim=c(0,0.15))
rb_sims$sims
mtext(side = 3, line = 0.25, adj = 0.5, paste("M =",M),cex = 1.5)
mtext(side = 3, line = -10, adj = 0.8, paste("RB Inefficiency =",eff_rb),cex = 1,col='red')
mtext(side = 3, line = -11, adj = 0.8, paste("Frac Inefficiency =",eff_frac),cex = 1,col='blue')
legend(x=0.3,y=15,c("RB","Fraction"),c("red","blue"))

quants = seq(0,0.95,by=0.1)
ymin = min(quantile(rb_sims$times,min(quants)),quantile(frac_sims$times,min(quants)))
ymax = max(quantile(rb_sims$times,max(quants)),quantile(frac_sims$times,max(quants)))

plot(quantile(rb_sims$times,quants),quants, col='red',type='l',xlim=c(ymin,ymax),ylab='Quantiles (0.05-0.95)',xlab='Time taken',
cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
lines(quantile(frac_sims$times,quants),quants,type='l',col='blue')
abline(v=mean(rb_sims$times),lty='dashed',col='red')
abline(v=mean(frac_sims$times),lty='dashed',col='blue')

dev.off()

### M = 100
pdf("simulation-study/BDI/N=100/different-m/M=100.pdf")
N = 1000
M = 100
frac_sims = list(sims=rep(0,N),times=rep(0,N))
rb_sims = list(sims=rep(0,N),times=rep(0,N))

for(i in 1:N){
  start.time = Sys.time()
  frac_sims$sims[i] = frac_list("BDI",x, theta, tout_list, obs_list,M)
  end.time = Sys.time()
  frac_sims$times[i] = end.time-start.time
}

for(i in 1:N){
  start.time = Sys.time()
  rb_sims$sims[i] = RB_list("BDI",x, theta, tout_list, lower_list, upper_list, tau,obs_list,M)
  end.time = Sys.time()
  rb_sims$times[i] = end.time-start.time
}
eff_rb = formatC(var(rb_sims$sims)*sum(rb_sims$times),format='e')
eff_frac = formatC(var(frac_sims$sims)*sum(frac_sims$times),format='e')

p1 = hist(frac_sims$sims,breaks=15,prob=T,col=rgb(0,0,1,1/4), xlim=c(0,0.15),ylim=c(0,40),
main='Estimates of likelihood',cex.lab=1.5,cex.main=1.5,cex.axis=1.5,xlab='likelihood',ylab='density')
p2 = hist(rb_sims$sims,breaks=15,freq=F,prob=T,add=T, col=rgb(1,0,0,1/4), xlim=c(0,0.15))
mtext(side = 3, line = 0.25, adj = 0.5, paste("M =",M),cex = 1.5)
mtext(side = 3, line = -10, adj = 0.8, paste("RB Inefficiency =",eff_rb),cex = 1,col='red')
mtext(side = 3, line = -11, adj = 0.8, paste("Frac Inefficiency =",eff_frac),cex = 1,col='blue')
legend(x=0.1,y=15,c("RB","Fraction"),c("red","blue"))

quants = seq(0,0.95,by=0.1)
ymin = min(quantile(rb_sims$times,min(quants)),quantile(frac_sims$times,min(quants)))
ymax = max(quantile(rb_sims$times,max(quants)),quantile(frac_sims$times,max(quants)))

plot(quantile(rb_sims$times,quants),quants, col='red',type='l',xlim=c(ymin,ymax),ylab='Quantiles (0.05-0.95)',xlab='Time taken',
cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
lines(quantile(frac_sims$times,quants),quants,type='l',col='blue')
abline(v=mean(rb_sims$times),lty='dashed',col='red')
abline(v=mean(frac_sims$times),lty='dashed',col='blue')

dev.off()

### M = 1000
pdf("simulation-study/BDI/N=100/different-m/M=1000.pdf")
N = 1000
M = 1000
frac_sims = list(sims=rep(0,N),times=rep(0,N))
rb_sims = list(sims=rep(0,N),times=rep(0,N))

for(i in 1:N){
  start.time = Sys.time()
  frac_sims$sims[i] = frac_list("BDI",x, theta, tout_list, obs_list,M)
  end.time = Sys.time()
  frac_sims$times[i] = end.time-start.time
}

for(i in 1:N){
  start.time = Sys.time()
  rb_sims$sims[i] = RB_list("BDI",x, theta, tout_list, lower_list, upper_list, tau,obs_list,M)
  end.time = Sys.time()
  rb_sims$times[i] = end.time-start.time
}
eff_rb = formatC(var(rb_sims$sims)*sum(rb_sims$times),format='e')
eff_frac = formatC(var(frac_sims$sims)*sum(frac_sims$times),format='e')

p1 = hist(frac_sims$sims,breaks=15,prob=T,col=rgb(0,0,1,1/4), xlim=c(0,0.15),ylim=c(0,100),
main='Estimates of likelihood',cex.lab=1.5,cex.main=1.5,cex.axis=1.5,xlab='likelihood',ylab='density')
p2 = hist(rb_sims$sims,breaks=15,freq=F,prob=T,add=T, col=rgb(1,0,0,1/4), xlim=c(0,0.15))
mtext(side = 3, line = 0.25, adj = 0.5, paste("M =",M),cex = 1.5)
mtext(side = 3, line = -10, adj = 0.8, paste("RB Inefficiency =",eff_rb),cex = 1,col='red')
mtext(side = 3, line = -11, adj = 0.8, paste("Frac Inefficiency =",eff_frac),cex = 1,col='blue')
legend(x=0.1,y=15,c("RB","Fraction"),c("red","blue"))

quants = seq(0,0.95,by=0.1)
ymin = min(quantile(rb_sims$times,min(quants)),quantile(frac_sims$times,min(quants)))
ymax = max(quantile(rb_sims$times,max(quants)),quantile(frac_sims$times,max(quants)))

plot(quantile(rb_sims$times,quants),quants, col='red',type='l',xlim=c(ymin,ymax),ylab='Quantiles (0.05-0.95)',xlab='Time taken',
cex.lab=1.5,cex.main=1.5,cex.axis=1.5)
lines(quantile(frac_sims$times,quants),quants,type='l',col='blue')
abline(v=mean(rb_sims$times),lty='dashed',col='red')
abline(v=mean(frac_sims$times),lty='dashed',col='blue')

dev.off()

