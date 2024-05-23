source("examples/Parameters/BDI/Drift/Drift (outlier -2SD).R")

### REMINDER: Reusing estimates?

M = 100
N = 10000
goal = 0.95

lower_limit_list = get_lower("BDI", obs_list,120)
upper_limit_list = get_upper("BDI", obs_list,120)

### tau = 0.10
pdf("simulation-study/BDI/N=100/special-case/outlier/2SD/different-drift/negative/different-tau/tau=0_10.pdf")
tau = 0.10
subtitle = paste('tau =',tau)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.0004,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### tau = 0.30
pdf("simulation-study/BDI/N=100/special-case/outlier/2SD/different-drift/negative/different-tau/tau=0_30.pdf")
tau = 0.30
subtitle = paste('tau =',tau)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.0004,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()

### tau = 0.60
pdf("simulation-study/BDI/N=100/special-case/outlier/2SD/different-drift/negative/different-tau/tau=0_60.pdf")
tau = 0.60
subtitle = paste('tau =',tau)
box = get_box("BDI",theta,tout_list,tau,x,obs_list,goal)
lower_list = box$lower
upper_list = box$upper

frac_sims = get_study_frac("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)
rb_sims = get_study_rb("BDI",x,theta,tout_list,obs_list,M,N,lower_list,upper_list,tau)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(50,150))
plot_study_estimates_generic(rb_sims,frac_sims,0.0004,0.99,subtitle)
plot_tables("BDI",lower_list,upper_list,theta,tau,tout_list,x,obs_list,lower_limit_list,upper_limit_list,M,N,rb_sims,frac_sims)
plot_study_estimates_M_time(rb_sims,frac_sims,30)
dev.off()
