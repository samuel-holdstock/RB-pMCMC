source("examples/Parameters/SIR_N/Equilibrium (N=1000).R")
lower = c(0,0)
upper = c(5,80)
tau = 20

M = 10000
print(RB("SIR_N", x, theta, tout, lower, upper, tau, obs, M))

1-pf(3.364471,1,7)
1-pf(0.31577,1,7)

(52.36-50.10)/(3-2)/(50.10/(10-4))

2*(1-pnorm(1.73))

128+19/2+c(-1,1)*2*(19/2 /qnorm(1-0.05/2)*qnorm(1-0.01/2))

pbinom(0,size=8,prob=11/38)
sum(dbinom(x=0,size=c(8:100),prob=11/38))

1-pgeom(7,11/38)


1-(1-(1-11/38)^8)
