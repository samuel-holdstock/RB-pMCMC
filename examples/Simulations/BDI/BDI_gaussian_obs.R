lambda = 1/5
mu = 3/10
gamma = 100

x = c(100)
theta = c(lambda,mu,gamma/10)

tout_list = seq(0,0.3,by=0.1)[-1]
tout_list = 1:10
obs_list = get_obs("BDI", x, theta, tout_list)
lower_list = get_lower("BDI",obs_list,c(8))
upper_list = get_upper("BDI",obs_list,c(8))
tau_list = get_tau_list("BDI",x,theta,obs_list,tout_list,20)

# Noisy Obs
NOISE = 1
obs_list_noise = obs_list + rnorm(prod(dim(obs_list)),mean=0,sd=NOISE)

plot_RB_MJP("BDI",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,50,c('Population'),c('blue'),c(0,500))
RB_particle_filter("BDI",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,500)

M = 2000
boot1 = rep(0,M)
RB_particles = 200
for(i in 1:M){
  boot1[i] = RB_particle_filter("BDI",x,theta,tout_list, lower_list, upper_list,tau_list,obs_list_noise,RB_particles)
}
var(boot1)
var(exp(boot1))/mean(exp(boot1))^2

boot2 = rep(0,M)
frac_particles=400
for(i in 1:M){
  boot2[i] = frac_particle_filter("BDI",x,theta,tout_list,obs_list,frac_particles)
}
var(boot2)
var(exp(boot2))/mean(exp(boot2))^2

get_benchmark = function(num_particles,times=100){
  bench_results = list()
  for(counter in 1:length(num_particles)){
    bench_results[[counter]] = microbenchmark(
      A={RB_particle_filter("BDI",x,theta,tout_list, floor(obs_list_noise-5), ceiling(obs_list_noise+5),tau_list,obs_list_noise,num_particles[counter])},
      B={RB_particle_filter("BDI",x,theta,tout_list, floor(obs_list_noise-1), ceiling(obs_list_noise+1),tau_list,obs_list_noise,num_particles[counter])},
      C={frac_particle_filter("BDI",x,theta,tout_list,obs_list,num_particles[counter])},
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
  A={cexpm(tau_list[1,]*get_coffin_matrix("BDI",floor(obs_list_noise-5)[1,], ceiling(obs_list_noise+5)[1,],theta))},
  B={cexpm(tau_list[2,]*get_coffin_matrix("BDI",floor(obs_list_noise-5)[2,], ceiling(obs_list_noise+5)[2,],theta))},
  C={cexpm(tau_list[3,]*get_coffin_matrix("BDI",floor(obs_list_noise-5)[3,], ceiling(obs_list_noise+5)[3,],theta))},
  D={cexpm(tau_list[1,]*get_coffin_matrix("BDI",floor(obs_list_noise-1)[1,], ceiling(obs_list_noise+1)[1,],theta))},
  E={cexpm(tau_list[2,]*get_coffin_matrix("BDI",floor(obs_list_noise-1)[2,], ceiling(obs_list_noise+1)[2,],theta))},
  F={cexpm(tau_list[3,]*get_coffin_matrix("BDI",floor(obs_list_noise-1)[3,], ceiling(obs_list_noise+1)[3,],theta))},
  unit='ms',times=500
))

plot(num_particles[-1],bench[[2]][1,][-1],type='o',ylim=c(0,15),col='red')
lines(num_particles[-1],bench[[2]][2,][-1],col='darkgreen',type='o')
lines(num_particles[-1],bench[[2]][3,][-1],col='blue',type='o')

abline(h=sum(second_bench$mean[1:3]),lty='dashed')
lines(num_particles,sum(second_bench$mean[1:3])+bench[[2]][3,])
lines(num_particles,second_bench$mean[4:6]+bench[[2]][3,])

time_RB = Sys.time()
r1 = mcmc_RB("BDI",x,theta, obs_list, lower_list, upper_list, tout_list, 1000, RB_particles,diag(3))
time_RB = Sys.time()-time_RB

time_frac = Sys.time()
r2 = mcmc_frac("BDI",x,theta, obs_list, tout_list, 1000, frac_particles,diag(3))
time_frac = Sys.time()-time_frac

r1$acceptance_rate
r2$acceptance_rate

time_RB = Sys.time()
for(i in 1:1000){
  RB_particle_filter("BDI",x,theta+rnorm(3,sd=0.001),tout_list,lower_list, upper_list,tau_list,obs_list_noise,RB_particles)
}
time_RB = Sys.time()-time_RB

time_frac = Sys.time()
for(i in 1:1000){
  frac_particle_filter("BDI",x,theta+rnorm(3,sd=0.001),tout_list,obs_list_noise,frac_particles)
}
time_frac = Sys.time()-time_frac

tau_list = get_tau_list("BDI",x,theta,obs_list,tout_list,10)
B = 50000
d1 = RB_particle_filter_debug("BDI",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,B)
d2 = frac_particle_filter_debug("BDI",x,theta,tout_list,obs_list,B)
d3 = RB_particle_filter_debug("BDI",x,theta,tout_list,box$lower,box$upper,tau_list,obs_list,B)

sum(d1[[1]]$particle_weights)^2/sum(d1[[1]]$particle_weights^2)
sum(d2[[1]]$weights)^2/sum(d2[[1]]$weights^2)

sum(d1[[1]]$weights[1:(length(d1[[1]]$weights)-B)])/sum(d1[[1]]$weights)

box = get_box_brownian_fast_list("BDI",theta,tout_list,tau_list,x,obs_list,beta=2)
get_tau_list("BDI",x,theta,obs_list,tout_list,alpha=10)
