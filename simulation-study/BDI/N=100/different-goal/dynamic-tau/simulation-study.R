source("examples/Parameters/BDI/Equilibrium (N=100).R")

### REMINDER: Reusing estimates?

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[1,1] = 100
M = 100
N = 10000

lower_limit_list = get_lower("BDI", obs_list,80)
upper_limit_list = get_upper("BDI", obs_list,80)

### goal = 0.90
pdf("simulation-study/BDI/N=100/different-goal/dynamic-tau/goal=0_90.pdf")
goal = 0.90
tau = get_tau("BDI", x, theta, obs, tout, sqrt(goal))
subtitle = paste('goal =',goal)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,sqrt(goal))
lower_list = box$lower-1
upper_list = box$upper+1

get_PVR_big_box_tau("BDI", x, theta, obs, tout, tau)
1-get_variance_brownian_fast("BDI", theta, tout, tau, x, obs, box$lower, box$upper)
sqrt(0.90)

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

signif(3.9665e-05*M,3)
signif(5.1236e-04*M,3)

0.00098596/M
0.00089989/M

### goal = 0.95
pdf("simulation-study/BDI/N=100/different-goal/dynamic-tau/goal=0_95.pdf")
goal = 0.95
tau = get_tau("BDI", x, theta, obs, tout, sqrt(goal))
subtitle = paste('goal =',goal)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,sqrt(goal))
lower_list = box$lower-1
upper_list = box$upper+1

get_PVR_big_box_tau("BDI", x, theta, obs, tout, tau)
1-get_variance_brownian_fast("BDI", theta, tout, tau, x, obs, box$lower-1, box$upper+1)
sqrt(0.95)

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### goal = 0.99
pdf("simulation-study/BDI/N=100/different-goal/dynamic-tau/goal=0_99.pdf")
goal = 0.99
tau = get_tau("BDI", x, theta, obs, tout, sqrt(goal))
subtitle = paste('goal =',goal)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,sqrt(goal))
lower_list = box$lower-1
upper_list = box$upper+1

get_PVR_big_box_tau("BDI", x, theta, obs, tout, tau)
1-get_variance_brownian_fast("BDI", theta, tout, tau, x, obs, box$lower-1, box$upper+1)
sqrt(0.99)

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()


0.0010041/100
0.00091455/100

1.8185e-05*100
5.0377e-04*100
