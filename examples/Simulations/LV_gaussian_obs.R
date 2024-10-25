pred = 100
prey = 100

alpha = 0.5
beta = 0.0025
gamma = 0.3

x = c(prey,pred)
theta = c(alpha,beta,gamma)

tout_list = c(1,2,3,4)
obs_list = get_obs("LV", x, theta, tout_list)
lower_list = get_lower("LV",obs_list,c(16,16))
upper_list = get_upper("LV",obs_list,c(16,16))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,100)

# Noisy Obs
NOISE = 2
obs_list_noise = obs_list + rnorm(prod(dim(obs_list_noise)),mean=0,sd=NOISE)

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,10,c('Prey','Predator'),c('blue','red'),c(0,500))

RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,500)

