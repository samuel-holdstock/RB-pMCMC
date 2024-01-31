source("examples/Parameters/BDI/Equilibrium (N=1000).R")

lower = c(990)
upper = c(1010)
tau = 0.1

M = 1000
print(RB("BDI",x, theta, tout, lower, upper, tau,obs,M))
