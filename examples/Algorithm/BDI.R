source("examples/Parameters/BDI/Equilibrium (N=100).R")

lower = c(990)
upper = c(1010)
tau = 0.1

lower = 50
upper = 150
M = 100000
tau = 1
(RB("BDI",x, theta, tout, lower, upper, tau,obs,M))
Q = get_coffin_matrix("BDI",lower,upper,theta)
P = expm(Q*tau)[,2]
P = matrix(P,ncol=1)
P
estimate = sim_data("BDI",x, theta, tout, lower, upper, tau)
get_estimate(estimate,lower,upper,P,obs)
(estimate)

state_to_index(obs,lower,upper)
