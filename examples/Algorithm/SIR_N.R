source("examples/Parameters/SIR_N/Equilibrium (N=1000).R")
lower = c(0,0)
upper = c(5,80)
tau = 20

M = 1000
print(RB("SIR_N", x, theta, tout, lower, upper, tau, obs, M))

