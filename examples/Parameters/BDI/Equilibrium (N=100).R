P = 100
lambda = 1/5
mu = 3/10
gamma = 100
tout = 1
theta = c(lambda,mu,gamma/10)
x = c(P)

obs = c(100)

lower_limit = 60
upper_limit = 140
get_mu("BDI",x,theta)
get_var("BDI",x,theta)
