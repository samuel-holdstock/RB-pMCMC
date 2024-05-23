P = 100
lambda = 1/5
mu = 3/10
gamma = 100+20*10
tout = 1
theta = c(lambda,mu,gamma/10)
x = c(P)

100+get_mu("BDI",x,theta)-5*sqrt(get_var("BDI",x,theta))

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[1,1] = 75

