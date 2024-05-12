source("examples/Parameters/BDI/Equilibrium (outlier).R")

### REMINDER: Reusing estimates?

tau = 0.30
M = 100
N = 10000

obs_list

### goal = 0.90
pdf("simulation-study/BDI/N=100/add-outlier/different-goal/tau=0.30/goal=0_90.pdf")
goal = 0.90
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_goal(rb_sims,frac_sims,0.001,0.99,goal)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### goal = 0.95
pdf("simulation-study/BDI/N=100/add-outlier/different-goal/tau=0.30/goal=0_95.pdf")
goal = 0.95
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_goal(rb_sims,frac_sims,0.001,0.99,goal)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### goal = 0.99
pdf("simulation-study/BDI/N=100/add-outlier/different-goal/tau=0.30/goal=0_99.pdf")
goal = 0.99
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_goal(rb_sims,frac_sims,0.0006,0.99,goal)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

