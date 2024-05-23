source("examples/Parameters/BDI/Equilibrium (N=100).R")

### REMINDER: Reusing estimates?

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[1,1] = 100
tau = 0.30
N = 10000
goal = 0.90

lower_limit_list = get_lower("BDI", obs_list,80)
upper_limit_list = get_upper("BDI", obs_list,80)

### M = 10
pdf("simulation-study/BDI/N=100/different-m/tau=0.30/M=10.pdf")
M = 10
subtitle = paste('M =',M)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.004,0.95,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### M = 100
pdf("simulation-study/BDI/N=100/different-m/tau=0.30/M=100.pdf")
M = 100
subtitle = paste('M =',M)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.95,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### M = 1000
pdf("simulation-study/BDI/N=100/different-m/tau=0.30/M=1000.pdf")
M = 1000
subtitle = paste('M =',M)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.0006,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

