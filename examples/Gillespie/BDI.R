source("examples/Parameters/BDI/Equilibrium (N=1000).R")

lower = c(990)
upper = c(1010)
tau = 0.1
obs_index = state_to_index(obs,lower,upper)
data = sim_data("BDI",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(920,1080))
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
box = get_box("BDI",theta,tout_list,tau,x,obs_list,0.9)
lower_list = box$lower
upper_list = box$upper
plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Population"),c("blue"),c(80,120))

get_box_brownian_fast("BDI",theta,tout,tau,x,obs,0.9)

abline(a=x[1],b=get_mu("BDI",x,theta)[1])

