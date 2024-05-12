source("examples/Parameters/BDI/Equilibrium (N=1000).R")

### REMINDER: Reusing estimates?

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
tau = 0.1
box = get_box("BDI",theta,tout_list,tau,x,obs_list,0.9)
lower_list = box$lower
upper_list = box$upper

pdf("simulation-study/BDI/N=1000/different-m/plot.pdf")
plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(900,1100))
dev.off()

obs_list
x

### M = 10
pdf("simulation-study/BDI/N=1000/different-m/M=10.pdf")
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

plot_study_estimates_M(rb_sims,frac_sims,0.002,1)
plot_study_estimates_M_time(rb_sims,frac_sims,30)

dev.off()

### M = 100
pdf("simulation-study/BDI/N=1000/different-m/M=100.pdf")
N = 10000
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

plot_study_estimates_M(rb_sims,frac_sims,0.001,0.9)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### M = 1000
pdf("simulation-study/BDI/N=1000/different-m/M=1000.pdf")
N = 10000
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

plot_study_estimates_M(rb_sims,frac_sims,0.0005,1)
plot_study_estimates_M_time(rb_sims,frac_sims,30)

dev.off()


