P = 1000
lambda = 1/5
mu = 3/10
gamma = 100+30
tout = 1
theta = c(lambda,mu,gamma)
x = c(P)

obs = c(1000)

lower_limit = 920
upper_limit = 1080
get_mu("BDI",x,theta)
get_var("BDI",x,theta)
lambda*P+mu*P+gamma
