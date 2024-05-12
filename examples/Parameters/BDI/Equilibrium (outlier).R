P = 100
lambda = 1/5
mu = 3/10
gamma = 100
tout = 1
theta = c(lambda,mu,gamma/10)
x = c(P)

obs = c(116)
lower_limit = 60
upper_limit = 140
get_mu("BDI",x,theta)
2*sqrt(get_var("BDI",x,theta))

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[1,1] = 116

