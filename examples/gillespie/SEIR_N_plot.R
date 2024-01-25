N = 1000
S0 = 997
E0 = 0
I0 = 3
R0 = N-S0-I0
beta = 0.4
alpha = 0.25
gamma = 0.04
total_time = 125

x = c(S0,E0,I0)
theta = c(N,beta,alpha,gamma)

lower = c(0,0,0)
upper = c(5,5,80)

obs = c(0,0,20)
obs_index = state_to_index(obs,lower,upper)
tau = 20

data = sim_data("SEIR_N",x,theta, total_time, lower, upper, 0.1)$data
data_chain = cbind(data,N - apply(data[,-1],1,sum))
plot_chain(data_chain,c("Susceptible","Exposed","Infected","Removed"),c("blue","orange","darkgreen","red"),ylim=c(0,1000))
add_box(lower,upper,total_time,tau,c("blue","orange","darkgreen"))
add_obs(total_time,obs,c("blue","orange","darkgreen"))

M = 1000
RB("SIR_N",x,theta, total_time, lower, upper, tau,obs,M)

