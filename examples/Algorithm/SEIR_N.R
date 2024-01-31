source("examples/Parameters/SEIR_N/Equilibrium (N=1000).R")
lower = c(0,0,0)
upper = c(5,5,80)
tau = 20

M = 1000
print(RB("SEIR_N", x, theta, tout, lower, upper, tau, obs, M))


