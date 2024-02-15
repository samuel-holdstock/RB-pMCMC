source("examples/Parameters/BDI/Equilibrium (N=1000).R")


lower = c(990)
upper = c(1010)
tau = 0.1
obs_index = state_to_index(obs,lower,upper)

data = sim_data("BDI",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(950,1050))
add_box(lower,upper,tout,tau,"blue")
add_obs(tout,obs,"blue")


obs = 1060
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(900,1100))
add_obs(tout,obs,"blue")
add_box_goal("BDI",0.9,theta,x,obs,tout,0.2,"blue")
data = sim_data("BDI",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
add_chain(data_chain,c("Population"),c("blue"))
