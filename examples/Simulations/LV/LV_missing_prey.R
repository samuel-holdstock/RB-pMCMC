pred = 100
prey = 100

alpha = 0.5
beta = 0.0025
gamma = 0.3

x = c(prey,pred)
theta = c(alpha,beta,gamma)

tout_list = c(1,2,3)*0.1
obs_list = get_obs("LV", x, theta, tout_list)
lower_list = get_lower("LV",obs_list,c(8,8))
upper_list = get_upper("LV",obs_list,c(8,8))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,20)

# Noisy Obs
NOISE = 1
obs_list_noise = obs_list + rnorm(prod(dim(obs_list)),mean=0,sd=NOISE)

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,100,c('Prey','Predator'),c('blue','red'),c(50,200))

prod((upper_list-lower_list)[1,]+1)+1

d1 = RB_particle_filter_debug("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,10000)
d2 = frac_particle_filter_debug("LV",x,theta,tout_list,obs_list_noise,10000)
d3 = RB_particle_filter_debug("LV",x,theta,tout_list,predicted_lower,predicted_upper,tau_list,obs_list_noise,10000)
c(d1[[1]]$exp_time*10^(-6),d1[[1]]$entire_time*10^(-6))
c(d2[[1]]$entire_time*10^(-6))

head(d1[[1]]$weights)

for(i in 1:length(tout_list)){
  print(c(paste0("Obs ",i),
    sum(d1[[i]]$particle_weights)^2/sum(d1[[i]]$particle_weights^2),
    sum(d2[[i]]$weights)^2/sum(d2[[i]]$weights^2),
    sum(d3[[i]]$particle_weights)^2/sum(d3[[i]]$particle_weights^2)
  ))
}
box_details(d1,tout_list)
box_details(d3,tout_list)

time = 3
t = if(time==1){tout_list[time]-tau_list[time]}else{tout_list[time]-tout_list[time-1]-tau_list[time]}
plot(jitter(d1[[time]]$xt_data),col='blue',pch='x')
points(jitter(d1[[time]]$xttau_data),xlim=c(90,130),ylim=c(90,130),col='red',pch='o')
V = diag(2)*1
predicted_mean = get_predicted_mean_matrix("LV",t,x,obs_list,tout_list,theta,V)
predicted_lower = get_predicted_lower_matrix("LV",t,x,obs_list,tout_list,theta,c(8,8),V)
predicted_upper = get_predicted_upper_matrix("LV",t,x,obs_list,tout_list,theta,c(8,8),V)
plot_box(obs_list,time,lower_list,upper_list,'darkgreen')
plot_box(predicted_mean,time,predicted_lower,predicted_upper,'darkmagenta')

M = 500
boot1 = rep(0,M)
RB_particles = 550
for(i in 1:M){
  boot1[i] = RB_particle_filter("LV",x,theta,tout_list, lower_list, upper_list,tau_list,obs_list_noise,RB_particles)
}
var(boot1)
var(exp(boot1))/mean(exp(boot1))^2

boot2 = rep(0,M)
frac_particles = 4000
for(i in 1:M){
  boot2[i] = frac_particle_filter("LV",x,theta,tout_list,obs_list_noise,frac_particles)
}
var(boot2)
var(exp(boot2))/mean(exp(boot2))^2

boot3 = rep(0,M)
RB_particles = 550
for(i in 1:M){
  boot3[i] = RB_particle_filter("LV",x,theta,tout_list, predicted_lower, predicted_upper,tau_list,obs_list_noise,RB_particles)
}
var(boot3)
var(exp(boot3))/mean(exp(boot3))^2


time_RB = Sys.time()
r1 = mcmc_RB("LV",x,theta, obs_list_noise, lower_list,upper_list, tout_list, 10000, RB_particles,diag(3))
time_RB = Sys.time()-time_RB

time_frac = Sys.time()
r2 = mcmc_frac("LV",x,theta, obs_list, tout_list, 10000, frac_particles,diag(3))
time_frac = Sys.time()-time_frac

time_RB2 = Sys.time()
r3 = mcmc_RB("LV",x,theta, obs_list_noise, predicted_lower,predicted_upper, tout_list, 10000, RB_particles,diag(3))
time_RB2 = Sys.time()-time_RB2

r1$acceptance_rate # 0.5583
r2$acceptance_rate # 0.4779
r3$acceptance_rate # 0.5609

effectiveSize(r1$ll) # 1685.759
effectiveSize(r2$ll) # 1510.364
effectiveSize(r3$ll) # 2004.543
