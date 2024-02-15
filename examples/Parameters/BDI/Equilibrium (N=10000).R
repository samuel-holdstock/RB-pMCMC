P = 10000
lambda = 1/5
mu = 3/10
gamma = 100
tout = 1
theta = c(lambda,mu,gamma*10)
x = c(P)

obs = c(10000)


get_mu("BDI",x,theta)
get_var("BDI",x,theta)
