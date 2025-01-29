pred = 100
prey = 100

alpha = 0.5
beta = 0.0025
gamma = 0.3

x = c(prey,pred)
theta = c(alpha,beta,gamma)

tout_list = c(1,2,3)
obs_list = get_obs("LV", x, theta, tout_list)
lower_list = get_lower("LV",obs_list,c(2,2))
upper_list = get_upper("LV",obs_list,c(2,2))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,5)

Q = get_coffin_matrix("LV",lower_list[1,]+1,upper_list[1,]-1,theta)
cexpm(Q*0.01)[,10]
cexpm(t(Q)*0.01)[10,]

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,100,c('Prey','Predator'),c('blue','red'),c(0,250))

RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,100)
frac_particle_filter("LV",x,theta,tout_list,obs_list,100)

d = RB_particle_filter_debug("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,10)

d[[1]]$xttau_data[d[[1]]$hitS==1,]
d[[1]]$xttau_data[d[[1]]$inS==1,]

sum(d[[1]]$hitS==1)/length(d[[1]]$hitS)
sum(d[[1]]$inS==1)/sum(d[[1]]$hitS==1)
sum(d[[1]]$inS==1)/length(d[[1]]$hitS)
Q = get_coffin_matrix("LV",lower_list[1,],upper_list[1,],theta)
plot(cexpm(Q*tau_list[1,])[,nrow(Q)],type='l')


plot(d[[1]]$xt_data+rnorm(length(d[[1]]$xt_data),sd=0.01))
points(t(obs_list[1,]),col='red',pch=4,cex=3)
lines(c(lower_list[1,1],lower_list[1,1]),c(lower_list[1,2],upper_list[1,2]),col='blue')
lines(c(lower_list[1,1],upper_list[1,1]),c(lower_list[1,2],lower_list[1,2]),col='blue')
lines(c(upper_list[1,1],upper_list[1,1]),c(lower_list[1,2],upper_list[1,2]),col='blue')
lines(c(lower_list[1,1],upper_list[1,1]),c(upper_list[1,2],upper_list[1,2]),col='blue')

boot1 = replicate(1000,RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,600))
boot2 = replicate(1000,frac_particle_filter("LV",x,theta,tout_list,obs_list,600))
30000*(1-(1-(1-exp(log(4226/30000)/9)))^(25))
(1-(1-exp(log(4226/30000)/9)))^(9)

hist(boot2)


microbenchmark(
  A={RB_particle_filter2("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,300)},
  B={frac_particle_filter2("LV",x,theta,tout_list,obs_list,400)},
  C={get_tau_list("LV",x,theta,obs_list,tout_list,4)},
  unit="ms",times=1000
)


var(exp(boot1))/mean(exp(boot1))^2
var(exp(boot2))/mean(exp(boot2))^2
hist(boot1)


sum(is.finite(boot1))
var(boot1[is.finite(boot1)])
sum(is.finite(boot2))
var(boot2[is.finite(boot2)])

time_RB = Sys.time()
r1 = mcmc_RB("LV",x,theta, obs_list, floor(obs_list-1), ceiling(obs_list+1), tout_list, 100, 150,diag(3))
# r1 = mcmc_RB("LV",x,theta, obs_list, lower_list+2, upper_list-2, tout_list, 1000, 100,diag(3))
time_RB = Sys.time()-time_RB

# c1 = cov(r1$log_posterior)
time_frac = Sys.time()
r2 = mcmc_frac("LV",x,theta, obs_list, tout_list, 100, 100,diag(3))
time_frac = Sys.time()-time_frac

r1$acceptance_rate*1000
r2$acceptance_rate*1000

plot(r1$log_posterior[,1],r1$log_posterior[,2],col='red')
points(r2$log_posterior[,1],r2$log_posterior[,2],col='blue')
points(log(theta[1]),log(theta[2]),pch=4,cex=4,lwd=3)

effectiveSize(r1$log_posterior)
effectiveSize(r2$log_posterior)

get_benchmark = function(num_particles,times=100){
  bench_results = list()
  for(counter in 1:length(num_particles)){
    bench_results[[counter]] = microbenchmark(
      A={RB_particle_filter("LV",x,theta,tout_list, floor(obs_list-5), ceiling(obs_list+5),tau_list,obs_list,num_particles[counter])},
      B={RB_particle_filter("LV",x,theta,tout_list, floor(obs_list-1), ceiling(obs_list+1),tau_list,obs_list,num_particles[counter])},
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
bench = get_benchmark(num_particles,100)
second_bench = summary(microbenchmark(
  A={cexpm(get_coffin_matrix("LV",floor(obs_list-5)[1,], ceiling(obs_list_noise+5)[1,],theta))},
  B={cexpm(get_coffin_matrix("LV",floor(obs_list-5)[2,], ceiling(obs_list_noise+5)[2,],theta))},
  C={cexpm(get_coffin_matrix("LV",floor(obs_list-5)[3,], ceiling(obs_list_noise+5)[3,],theta))},
  E={index_to_state(12,floor(obs_list-5)[1,], ceiling(obs_list_noise+5)[1,])},
  F={state_to_index(obs_list[1,],floor(obs_list-5)[1,], ceiling(obs_list+5)[1,])},
  H={cexpm(get_coffin_matrix("LV",floor(obs_list-1)[1,], ceiling(obs_list+1)[1,],theta))},
  I={cexpm(get_coffin_matrix("LV",floor(obs_list-1)[2,], ceiling(obs_list+1)[2,],theta))},
  J={cexpm(get_coffin_matrix("LV",floor(obs_list-1)[3,], ceiling(obs_list+1)[3,],theta))},
  unit='ms',times=1000
))

plot(num_particles[-1],bench[[2]][1,][-1],type='o',ylim=c(0,15),col='red')
lines(num_particles[-1],bench[[2]][2,][-1],col='darkgreen',type='o')
lines(num_particles[-1],bench[[2]][3,][-1],col='blue',type='o')
lines(num_particles[-1],bench[[2]][4,][-1],col='purple',type='o')
lines(num_particles[-1],bench[[2]][5,][-1],col='orange',type='o')

sum(d[[1]]$weights[1:9])
sum(d[[1]]$weights)

tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,30)
d1 = RB_particle_filter_debug("LV",x,theta,tout_list,lower_list-5,upper_list+5,tau_list,obs_list,10000)
d2 = frac_particle_filter_debug("LV",x,theta,tout_list,obs_list,10000)
Q = get_coffin_matrix("LV",lower_list[1,],upper_list[1,],theta)
cexpm(Q*tau_list[1,])

sum(d1[[1]]$particle_weights)^2/sum(d1[[1]]$particle_weights^2)
sum(d2[[1]]$weights)^2/sum(d2[[1]]$weights^2)

d1[[1]]$presample[order(d1[[1]]$weights),]
d1[[1]]$weights[order(d1[[1]]$weights),]

mean(d1[[1]]$inS)

d1[[1]]$presample[which.min(d1[[1]]$weights),]

d2[[1]]$presample[order(d2[[1]]$weights),]

obs_list[1,]

get_f(c(132,94),obs_list[1,])


sum(which(d1[[1]]$weights<1e-20))
sum(which(d1[[2]]$weights<1e-20))
