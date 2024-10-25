source("examples/Parameters/BDI/Equilibrium (N=100).R")
source("examples/Parameters/BDI/Drift/Drift (N=100).R")
get_mu("BDI",x,theta)
get_covar("BDI",x,theta)

microbenchmark(
  A={state_to_index(x,lower,upper)},
  B={index_to_state(10,lower,upper)},
  C={RB("BDI",x,theta,tout,lower,upper,tau,obs,1)},
  times=10000
)

get_probability_hit(x,obs,tout,20,sqrt(get_covar("BDI",x,theta)))
get_density_hit(x,obs,tout,20,sqrt(get_covar("BDI",x,theta)))
get_MVT_density_hit(x,obs,tout,get_mu("BDI",x,theta),get_covar("BDI",x,theta))
get_VRF_big_box_tau("BDI",x,theta,obs,tout,0.9999)
get_VRF_big_box_tau2("BDI",x,theta,obs,tout,0.9999)

get_tau("BDI",x,theta,obs,tout,30)
get_tau2("BDI",x,theta,obs,tout,30)

p = get_MVT_density_hit(x,obs,tout,get_mu("BDI",x,theta),get_covar("BDI",x,theta))
E = monte_carlo_var_inf("BDI", x, theta, obs, tout, tau, 0.0001, 10000)
p*(1-p)/var(E)


Q = t(get_coffin_matrix("BDI",lower,upper,theta))
obs_index = state_to_index(obs,lower,upper)
v = rep(0,nrow(Q))
v[obs_index] = 1
vT_exp_Q(v,Q*tau,1e-20,F,T,F)
SS_exp_Q(Q*tau,1e-10,F)


c(t(P))
(expm(Q*tau)[9,]+expm(Q*tau)[10,]+expm(Q*tau)[11,])/3

sum(vT_exp_Q(v,Q*tau,1e-20,F,T,F))

expm(Q)[,22]

lower=60
upper=140
microbenchmark(
  A={get_coffin_matrix("BDI",lower,upper,theta)},
  B={index_to_state(10,lower,upper)},
  C={state_to_index(x,lower,upper)},
  D={vT_exp_Q(v,Q*tau,1e-20,F,T,F)},
  E={RB("BDI",x,theta,tout,lower,upper,tau,obs,1)}
)
  
lower = c(90)
upper = c(110)
tau = 0.1
obs_index = state_to_index(obs,lower,upper)
data = sim_data("BDI",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(50,150))
add_box(lower,upper,tout,tau,"blue")
add_obs(tout,obs,"blue")


obs = 1060
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(900,1100))
add_obs(tout,obs,"blue")
add_box_goal("BDI",0.9,theta,x,obs,tout,0.2,"blue")
data = sim_data("BDI",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
add_chain(data_chain,c("Population"),c("blue"))


theta = c(0.3,0.3)
obs = 1000
data = sim_data("BDI_2",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(920,1080))
add_box(lower,upper,tout,tau,"blue")
add_obs(tout,obs,"blue")

theta = c(0.2,0.3,70)
tout_list = c(1)
tout_list = c(0.25,0.5,0.75,1)
obs_list = get_obs("BDI", x, theta, tout_list)
# lower_list = get_lower("BDI",obs_list,c(10))
# upper_list = get_upper("BDI",obs_list,c(10))
box = get_box("BDI",theta,tout_list,tau,x,obs_list,100)
lower_list = box$lower
upper_list = box$upper

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,0.05,obs_list,30,c("Population"),c("blue"),c(50,150))
tau=0.25
tau_list

z = sim_data("BDI",x,theta, tout, 97, 103, 0.03)
z$xttau_data

z