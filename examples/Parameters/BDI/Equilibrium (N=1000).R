P = 1000
lambda = 1/5
mu = 3/10
gamma = 100
tout = 1
theta = c(lambda,mu,gamma)
x = c(P)

obs = c(1000)


get_mu("BDI",x,theta)
get_var("BDI",x,theta)