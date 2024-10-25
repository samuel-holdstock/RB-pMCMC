source("examples/Parameters/BDI/Equilibrium (N=100).R")

lower = c(97)
upper = c(103)
tau = 1

lower = 50
upper = 150
M = 100000
tau = 1
# 0.05417058

lower=50
upper = 150
tau = 1

RB("BDI",x,theta,tout,lower,upper,tau,obs,10000)
RB("BDI",x,theta,tout,lower,upper,0,obs,1000000)
RB_noise("BDI",x,theta,tout,95,105,0.3,obs,2)


tau_list = get_tau_list("BDI",x,theta,obs_list,tout_list,10)
lower_list = get_lower("BDI",obs_list,6)
upper_list = get_upper("BDI",obs_list,6)


RB_particle_filter("BDI", x, theta, tout_list,lower_list,upper_list, tau_list,obs_list,100000)
RB_particle_filter("BDI", x, theta, tout_list,lower_list,upper_list, tau_list*0,obs_list,100000)
RB_particle_filter("BDI", x, theta, tout_list,lower_list,upper_list, tout_list,obs_list,100000)


lower_list[1,1] = 97
upper_list[1,1] = 103
tau_list[1,1]=0.06


p_RB

p_frac = frac_noise("BDI",x,theta,tout,obs,500000)
p_frac

Q = get_coffin_matrix("BDI",lower,upper,theta)
P = expm(Q)
get_fP(lower,upper,obs,P)


sim_data("BDI",100,theta,tout,130,150,1)$xttau_data
RB("BDI",100,theta,tout,lower,upper,1,obs,100)
RB("BDI",100,theta,tout,99,150,1,obs,100000)

RB_particle_filter("BDI",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,100000)
RB_particle_filter("BDI",x,theta,tout_list,lower_list,upper_list,tau_list*0,obs_list,100000)
RB_particle_filter("BDI",x,theta,tout_list,lower_list,upper_list,rep(0.25,4),obs_list,100000)

