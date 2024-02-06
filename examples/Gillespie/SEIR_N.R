source("examples/Parameters/SEIR_N/Equilibrium (N=1000).R")
lower = c(1,0,0)
upper = c(5,5,80)

obs_index = state_to_index(obs,lower,upper)
tau = 20

data = sim_data("SEIR_N",x,theta, tout, lower, upper, tau)$data
data_chain = cbind(data,N - apply(data[,-1],1,sum))
plot_chain(data_chain,c("Susceptible","Exposed","Infected","Removed"),c("blue","orange","darkgreen","red"),ylim=c(0,1000))
add_box(lower,upper,tout,tau,c("blue","orange","darkgreen"))
add_obs(tout,obs,c("blue","orange","darkgreen"))

