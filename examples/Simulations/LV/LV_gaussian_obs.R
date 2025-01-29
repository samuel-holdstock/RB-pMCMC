pred = 100
prey = 100

alpha = 0.5
beta = 0.0025
gamma = 0.3

x = c(prey,pred)
theta = c(alpha,beta,gamma)

tout_list = seq(0,0.3,by=0.1)[-1]
tout_list = c(1,2,3)
obs_list = get_obs("LV", x, theta, tout_list)
lower_list = get_lower("LV",obs_list,c(8,8))
upper_list = get_upper("LV",obs_list,c(8,8))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,5)

# Noisy Obs
NOISE = 1
obs_list_noise = obs_list + rnorm(prod(dim(obs_list)),mean=0,sd=NOISE)

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,50,c('Prey','Predator'),c('blue','red'),c(0,500))
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,500)

M = 2000
boot1 = rep(0,M)
for(i in 1:M){
  boot1[i] = RB_particle_filter("LV",x,theta,tout_list, lower_list, upper_list,tau_list,obs_list_noise,100)
}
var(boot1)
var(exp(boot1))/mean(exp(boot1))^2

boot2 = rep(0,M)
for(i in 1:M){
  boot2[i] = frac_particle_filter("LV",x,theta,tout_list,obs_list,200)
}
var(boot2)
var(exp(boot2))/mean(exp(boot2))^2

get_benchmark = function(num_particles,times=100){
  bench_results = list()
  for(counter in 1:length(num_particles)){
    bench_results[[counter]] = microbenchmark(
      A={RB_particle_filter("LV",x,theta,tout_list, floor(obs_list_noise-5), ceiling(obs_list_noise+5),tau_list,obs_list_noise,num_particles[counter])},
      B={RB_particle_filter("LV",x,theta,tout_list, floor(obs_list_noise-1), ceiling(obs_list_noise+1),tau_list,obs_list_noise,num_particles[counter])},
      C={frac_particle_filter("LV",x,theta,tout_list,obs_list,num_particles[counter])},
      unit='ms',times=times
    )
  }
  bench_lower = sapply(bench_results,function(x){summary(x)$lq})
  bench_mean = sapply(bench_results,function(x){summary(x)$median})
  bench_upper = sapply(bench_results,function(x){summary(x)$uq})
  return(list(bench_lower,bench_mean,bench_upper))
}

num_particles = c(seq(5,25,by=5),seq(30,100,by=10))
bench = get_benchmark(num_particles,500)
second_bench = summary(microbenchmark(
  A={cexpm(tau_list[1,]*get_coffin_matrix("LV",floor(obs_list_noise-5)[1,], ceiling(obs_list_noise+5)[1,],theta))},
  B={cexpm(tau_list[2,]*get_coffin_matrix("LV",floor(obs_list_noise-5)[2,], ceiling(obs_list_noise+5)[2,],theta))},
  C={cexpm(tau_list[3,]*get_coffin_matrix("LV",floor(obs_list_noise-5)[3,], ceiling(obs_list_noise+5)[3,],theta))},
  D={cexpm(tau_list[1,]*get_coffin_matrix("LV",floor(obs_list_noise-1)[1,], ceiling(obs_list_noise+1)[1,],theta))},
  E={cexpm(tau_list[2,]*get_coffin_matrix("LV",floor(obs_list_noise-1)[2,], ceiling(obs_list_noise+1)[2,],theta))},
  F={cexpm(tau_list[3,]*get_coffin_matrix("LV",floor(obs_list_noise-1)[3,], ceiling(obs_list_noise+1)[3,],theta))},
  unit='ms',times=500
))

plot(num_particles[-1],bench[[2]][1,][-1],type='o',ylim=c(0,15),col='red')
lines(num_particles[-1],bench[[2]][2,][-1],col='darkgreen',type='o')
lines(num_particles[-1],bench[[2]][3,][-1],col='blue',type='o')

abline(h=sum(second_bench$mean[1:3]),lty='dashed')
lines(num_particles,sum(second_bench$mean[1:3])+bench[[2]][3,])
lines(num_particles,second_bench$mean[4:6]+bench[[2]][3,])

time_RB = Sys.time()
r1 = mcmc_RB("LV",x,theta, obs_list, lower_list, upper_list, tout_list, 1000, 40,diag(3))
time_RB = Sys.time()-time_RB

time_frac = Sys.time()
r2 = mcmc_frac("LV",x,theta, obs_list, tout_list, 1000, 100,diag(3))
time_frac = Sys.time()-time_frac

r1$acceptance_rate
r2$acceptance_rate

time_RB = Sys.time()
for(i in 1:1000){
  RB_particle_filter("LV",x,theta+rnorm(3,sd=0.001),tout_list,lower_list, upper_list,tau_list,obs_list_noise,100)
}
time_RB = Sys.time()-time_RB

time_frac = Sys.time()
for(i in 1:1000){
  frac_particle_filter("LV",x,theta+rnorm(3,sd=0.001),tout_list,obs_list_noise,200)
}
time_frac = Sys.time()-time_frac

