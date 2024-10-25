source("examples/Parameters/SIR_N/Equilibrium (N=1000).R")
lower = c(0,0)
upper = c(5,80)

obs_index = state_to_index(obs,lower,upper)
tau = 20

data = sim_data("SIR_N",x, theta, tout, lower, upper, 0.1)$data
data = cbind(data,N - apply(data[,-1],1,sum))
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Susceptible","Infected","Removed"),c("blue","darkgreen","red"),ylim=c(0,1000))
add_box(lower,upper,tout,tau,c("blue","darkgreen"))
add_obs(tout,obs,c("blue","darkgreen"))
