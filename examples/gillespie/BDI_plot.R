P = 1000
lambda = 0.1
gamma = 0.1
mu = 3
total_time = 1
theta = c(lambda,gamma,mu)
x = c(P)

lower = c(985)
upper = c(1015)

obs = c(1000)
obs_index = state_to_index(obs,lower,upper)
tau = 0.3

data = sim_data("BDI",x,theta, total_time, lower, upper, 0.1)$data
data_chain = get_chain(data,total_time)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(975,1050))
add_box(lower,upper,total_time,tau,"blue")
add_obs(total_time,obs,"blue")

M = 1000
RB("BDI",x, theta, total_time, lower, upper, tau,obs,M)

