pred = 200
prey = 200

alpha = 0.5
beta = 0.0025
gamma = 0.3
tout = 20

#alpha/beta
#gamma/beta

x = c(pred,prey)
theta = c(alpha,beta,gamma)

obs = c(292,112)
lower = c(250,70)
upper = c(330,150)

lower_limit = c(0,0)
upper_limit = c(30,30)

tout = 1
obs = c(220,190)

tout_list = c(1)
obs_list = get_obs("LV", x, theta, tout_list)
lower = c(obs[1]-30,obs[2]-30)
upper = c(obs[1]+30,obs[2]+30)

lower_list = get_lower("LV",obs_list,c(16,16))
upper_list = get_upper("LV",obs_list,c(16,16))
tau_list = get_tau("LV",x,theta,obs_list,tout_list,10)
tau_list = matrix(c(tau_list))
tau_list=matrix(c(0.06,0.06))
