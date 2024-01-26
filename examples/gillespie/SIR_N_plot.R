N = 1000
S0 = 997
I0 = 3
R0 = N-S0-I0
beta = 0.4
gamma = 0.04
total_time = 125

x = c(S0,I0)
theta = c(N,beta,gamma)

lower = c(0,0)
upper = c(5,80)

obs = c(0,13)
obs_index = state_to_index(obs,lower,upper)
tau = 20

data = sim_data("SIR_N",x, theta, total_time, lower, upper, 0.1)$data
data = cbind(data,N - apply(data[,-1],1,sum))
data_chain = get_chain(data,total_time)
plot_chain(data_chain,c("Susceptible","Infected","Removed"),c("blue","darkgreen","red"),ylim=c(0,1000))
add_box(lower,upper,total_time,tau,c("blue","darkgreen"))
add_obs(total_time,obs,c("blue","darkgreen"))

M = 1000
print(RB("SIR_N",x,theta, total_time, lower, upper, tau,obs,M))


