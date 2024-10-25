source("examples/Parameters/LV/(N=100).R")

get_mu("LV",x,theta)
get_covar("LV",x,theta)
get_MVT_density_hit(x,obs,tout,get_mu("LV",x,theta),get_covar("LV",x,theta))
get_VRF_big_box_tau2("LV",x,theta,obs,tout,0.1)
get_tau2("LV",x,theta,obs,tout,100)

tau = get_tau2("LV",x,theta,obs,tout,1000)
tau = 0.1
p = get_MVT_density_hit(x,obs,tout,get_mu("LV",x,theta),get_covar("LV",x,theta))
E = monte_carlo_var_inf("LV", x, theta, obs, tout, tau, 0.0001, 10000)
p*(1-p)/var(E)


noise_lower = c(195,190)
x
noise_upper = c(205,210)
noise_lower_index = state_to_index(noise_lower,lower,upper)
state_to_index(x,lower,upper)
state_to_index(noise_upper,lower,upper)



# 335 - (335+(220-190))
index_to_state(noise_lower_index,lower,upper)
index_to_state(noise_lower_index+(noise_upper[2]-noise_lower[2]),lower,upper)

index_to_state(noise_lower_index+(upper[2]-noise_lower[2]+1)+(noise_lower[2]-lower[2]),lower,upper)
index_to_state(noise_lower_index+(upper[2]-noise_lower[2]+1)+(noise_lower[2]-lower[2])+(noise_upper[2]-noise_lower[2]),lower,upper)

index_to_state(noise_lower_index+2*(upper[2]-noise_lower[2]+1)+2*(noise_lower[2]-lower[2]),lower,upper)
index_to_state(noise_lower_index+2*(upper[2]-noise_lower[2]+1)+2*(noise_lower[2]-lower[2])+(noise_upper[2]-noise_lower[2]),lower,upper)


mvrnormArma2(10,C)

plot(A,col='blue')
points(B,col='green')
C[1,1]=100

theta

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
box = get_box("LV",theta,tout_list,tau,x,obs_list,5)
lower_list = box$lower
upper_list = box$upper

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau,obs_list,30,c("Predator","Prey"),c("red","blue"),c(50,400))
tau = 0.06
as.numeric((upper_list - obs_list)[,1])
as.numeric((obs_list - lower_list)[,1])

abline(a=x[1],b=get_mu("LV",x,theta)[1])
abline(a=x[2],b=get_mu("LV",x,theta)[2])

tout_list = c(1,2,3,4)
obs_list = get_obs("LV", x, theta,tout_list)
obs_list_noise = obs_list+rnorm(length(obs_list),0,2)

lower_list = get_lower("LV",obs_list,rep(10,nrow(obs_list)))
upper_list = get_upper("LV",obs_list,rep(10,nrow(obs_list)))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,100)
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,c(tau_list),obs_list,100000)
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,c(tau_list*0),obs_list,100000)
boot = replicate(500,RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,c(tau_list),obs_list_noise,500))

lower_list = get_lower("LV",obs_list,c(1,1,1,1))
upper_list = get_upper("LV",obs_list,c(1,1,1,1))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,100)*1
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,c(tau_list),obs_list,100)
boot2 = replicate(500,RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,c(tau_list*0),obs_list_noise,500))

RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,c(tau_list),obs_list_noise,10000)
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,rep(1,nrow(obs_list)),obs_list_noise,1000000)
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,rep(0,nrow(obs_list)),obs_list_noise,1000000)
BSfilter(50,y,simX0,simNextX,logLik,thetas)$lml
1
tout_list=c(5,10)

y
obs_list_noise

var(boot)
hist(boot)
var(boot2)
hist(boot2)

mean(boot2)
mean(boot)




c = sim_data("LV",x,thetas,1,lower_list[1,],upper_list[1,],0)
c$inS
c$xt_data

lower_list[1,]
upper_list[1,]
