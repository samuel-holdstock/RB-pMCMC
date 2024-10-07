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
get_covar("BDI",x,theta)

gamma*(1+2*mu)/((1-(lambda-mu)^2)*(1-(lambda-mu)))

tout_list = c(1)
obs_list = get_obs("BDI", x, theta, tout_list)
obs_list[1,1] = 100
tau = 0.10

lower_limit_list = get_lower("BDI", obs_list,80)
upper_limit_list = get_upper("BDI", obs_list,80)
