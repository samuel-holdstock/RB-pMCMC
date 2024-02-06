source("examples/Parameters/SIR_N/Equilibrium (N=1000).R")

lower = c(0,0)
upper = c(5,80)
tau = 20

get_variance_brownian("SIR_N",theta,tout,tau,x,obs,lower,upper)
get_variance_brownian_fast("SIR_N",theta,tout,tau,x,obs,lower,upper)

index_to_state(,lower,upper)

x + tout*get_mu("SIR_N",x,theta)
get_var("SIR_N",x,theta)

