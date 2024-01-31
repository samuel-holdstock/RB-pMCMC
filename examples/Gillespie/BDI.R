source("examples/Parameters/BDI/Equilibrium (N=1000).R")

lower = c(990)
upper = c(1010)
tau = 0.1
obs_index = state_to_index(obs,lower,upper)

data = sim_data("BDI",x,theta, tout, lower, upper, tau)$data
data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(975,1050))
add_box(lower,upper,tout,tau,"blue")
add_obs(tout,obs,"blue")


get_variance_brownian("BDI", theta, tout, 0.1, x, obs, 1000,1000)

