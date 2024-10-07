source("examples/Parameters/LV/(N=100).R")

get_mu("LV",x,theta)
get_covar("LV",x,theta)
get_MVT_density_hit(x,obs,tout,get_mu("LV",x,theta),get_covar("LV",x,theta))
get_VRF_big_box_tau2("LV",x,theta,obs,tout,0.1)
get_tau2("LV",x,theta,obs,tout,100)

tau = get_tau2("LV",x,theta,obs,tout,10)
tau = 0.1
p = get_MVT_density_hit(x,obs,tout,get_mu("LV",x,theta),get_covar("LV",x,theta))
E = monte_carlo_var_inf("LV", x, theta, obs, tout, tau, 0.0001, 10000)
p*(1-p)/var(E)

obs_index = state_to_index(obs,lower,upper)
data = sim_data("LV",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Predator","Prey"),c("red","blue"),ylim=c(0,400))
add_box(lower,upper,tout,tau,c("red","blue"))
add_obs(tout,obs,c("red","blue"))

data = sim_data("LV",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
add_chain(data_chain,c("Predator","Prey"),c("red","blue"),ylim=c(0,400))

tau = 0.1
tout_list = c(1,2,3,4)
obs_list = get_obs("LV", x, theta,tout_list)
# lower_list = get_lower("LV",obs_list,c(10,10))
# upper_list = get_upper("LV",obs_list,c(10,10))
box = get_box("LV",theta,tout_list,tau,x,obs_list,0.95)
lower_list = box$lower
upper_list = box$upper
get_box_brownian_fast("LV",theta,tout,tau,x,obs,0.9)

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Predator","Prey"),c("red","blue"),c(50,400))

as.numeric((upper_list - obs_list)[,1])
as.numeric((obs_list - lower_list)[,1])

abline(a=x[1],b=get_mu("LV",x,theta)[1])
abline(a=x[2],b=get_mu("LV",x,theta)[2])
