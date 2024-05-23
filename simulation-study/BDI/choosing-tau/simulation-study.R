### PVR = 0.90
pdf("simulation-study/BDI/choosing-tau/goal=90.pdf")
PVR = 0.90
M = 100
N = 10000

source('examples/Parameters/BDI/Equilibrium (N=100).R')
plot_var_big_box_tau("BDI",x,theta,obs_list[1,1],tout_list[1])
plot_PVR_big_box_tau("BDI",x,theta,obs_list[1,1],tout_list[1])
plot_PVR_target_goal_brownian("BDI",x,theta,obs_list[1,1],tout_list[1])

goal = sqrt(PVR)
tau = get_tau("BDI",x,theta,obs_list[1,1],tout_list[1],goal)
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

### outlier, goal = 0.90
pdf("simulation-study/BDI/choosing-tau/outlier-goal=90.pdf")
PVR = 0.90
M = 100
N = 10000

source("examples/Parameters/BDI/Equilibrium (outlier 2SD).R")
plot_var_big_box_tau("BDI",x,theta,obs_list[1,1],tout_list[1])
plot_PVR_big_box_tau("BDI",x,theta,obs_list[1,1],tout_list[1])
plot_PVR_target_goal_brownian("BDI",x,theta,obs_list[1,1],tout_list[1])

goal = sqrt(PVR)
tau = get_tau("BDI",x,theta,obs_list[1,1],tout_list[1],goal)
subtitle = paste('tau =',tau)
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

