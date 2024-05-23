source("examples/Parameters/BDI/Equilibrium (N=100).R")

### REMINDER: Reusing estimates?

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[1,1] = 100
tau = 0.60
M = 100
N = 10000

lower_limit_list = get_lower("BDI", obs_list,80)
upper_limit_list = get_upper("BDI", obs_list,80)

### goal = 0.90
pdf("simulation-study/BDI/N=100/different-goal/tau=0.60/goal=0_90.pdf")
goal = 0.90
subtitle = paste('goal =',goal)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### goal = 0.95
pdf("simulation-study/BDI/N=100/different-goal/tau=0.60/goal=0_95.pdf")
goal = 0.95
subtitle = paste('goal =',goal)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### goal = 0.99
pdf("simulation-study/BDI/N=100/different-goal/tau=0.60/goal=0_99.pdf")
goal = 0.99
subtitle = paste('goal =',goal)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.002,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

