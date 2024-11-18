pred = 100
prey = 100

alpha = 0.5
beta = 0.0025
gamma = 0.3

x = c(prey,pred)
theta = c(alpha,beta,gamma)

tout_list = c(1,2,3)
obs_list = get_obs("LV", x, theta, tout_list)
lower_list = get_lower("LV",obs_list,c(16,16))
upper_list = get_upper("LV",obs_list,c(16,16))
tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,100)

# Noisy Obs
NOISE = 2
obs_list_noise = obs_list + rnorm(prod(dim(obs_list)),mean=0,sd=NOISE)

plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,50,c('Prey','Predator'),c('blue','red'),c(0,500))

RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,500)
RB_list("LV",x, theta, tout_list, lower_list, upper_list, tau_list, obs_list,500)

boot = replicate(100,RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,250))
var(boot)
plot(density(boot))

lower_list = get_lower("LV",obs_list,c(30,30))
upper_list = get_upper("LV",obs_list,c(30,30))

upper_list[2,2] = obs_list[2,2]+20
lower_list[2,2] = obs_list[2,2]-10
upper_list[3,2] = obs_list[3,2]+5
lower_list[3,2] = obs_list[3,2]-25
upper_list[4,2] = obs_list[4,2]+25
lower_list[4,2] = obs_list[4,2]-5

upper_list[,1]-lower_list[,1]
upper_list[,2]-lower_list[,2]

tau_list = get_tau_list("LV",x,theta,obs_list,tout_list,200)
plot_RB_MJP("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list,10,c('Prey','Predator'),c('blue','red'),c(0,500))

z3 = mcmc_RB("LV", x, theta, obs_list_noise, lower_list,upper_list, tout_list, 1000, 150, covar2)
covar = cov(z$log_posterior)

plot(density(exp(z$log_posterior[,1])))

theta
c(mean(exp(z$log_posterior[,1])),mean(exp(z$log_posterior[,2])),mean(exp(z$log_posterior[,3])))
plot(z$log_posterior[,1],type='l')

effectiveSize(z2$log_posterior)
effectiveSize(z3$log_posterior)
sqrt(2.43)

plot(z2$log_posterior[,1],type='l')

effectiveSize(apply(z2$log_posterior,MARGIN = 2,FUN = inv_transform))

M = 500
boot1 = rep(0,M)
for(i in 1:M){
  boot1[i] = RB_particle_filter("LV",x,theta,tout_list, floor(obs_list_noise-5), ceiling(obs_list_noise+5),tau_list,obs_list_noise,50)
}
var(boot1)
boot2 = rep(0,M)
for(i in 1:M){
  boot2[i] = RB_particle_filter("LV",x,theta,tout_list, floor(obs_list_noise-1), ceiling(obs_list_noise+1),tau_list,obs_list_noise,120)
}
var(boot2)

boot3 = rep(0,M)
for(i in 1:M){
  boot3[i] = frac_particle_filter("LV",x,theta,tout_list,obs_list,120)
}
var(boot3)

microbenchmark(
  A={RB_particle_filter("LV",x,theta,tout_list, floor(obs_list_noise-5), ceiling(obs_list_noise+5),tau_list,obs_list_noise,50)},
  B={RB_particle_filter("LV",x,theta,tout_list, floor(obs_list_noise-1), ceiling(obs_list_noise+1),tau_list*0,obs_list_noise,120)},
  C={frac_particle_filter("LV",x,theta,tout_list,obs_list,120)}
)

mean(boot1)
mean(boot2)
mean(boot3)+c(-2,2)*sd(boot3)

sd(boot1)
sd(boot2)
sd(boot3)

plot(density(boot1),col='blue',xlim=c(-30,-36))
lines(density(boot2),col='red')
lines(density(boot3),col='darkgreen')

plot_RB_MJP("LV",x,theta,tout_list,lower_list-30,upper_list+30,tau_list,obs_list,10,c('Prey','Predator'),c('blue','red'),c(0,500))


RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,50000)
RB_particle_filter("LV",x,theta,tout_list,lower_list,upper_list,tau_list*0,obs_list_noise,50000)

frac_particle_filter("LV",x,theta,tout_list,obs_list_noise,50000)

RB_particle_filter_gibbs("LV",x,theta,tout_list,lower_list,upper_list,tau_list,obs_list_noise,50)
frac_particle_filter_gibbs("LV",x,theta,tout_list,obs_list_noise,50)


z1 = mcmc_RB_marginal("LV",x,theta, obs_list_noise, floor(obs_list_noise-5), ceiling(obs_list_noise+5), tout_list, 1000, 100)
z2 = mcmc_RB_marginal("LV",x,theta, obs_list_noise, floor(obs_list_noise-1), ceiling(obs_list_noise+1), tout_list, 1000, 150)
z3 = mcmc_frac_marginal("LV",x,theta, obs_list_noise, tout_list, 1000, 150)

microbenchmark(expm(matrix(runif(9),nrow=3,ncol=3)),unit='s')

z1$acceptance_rate
z2$acceptance_rate
z3$acceptance_rate


mcmc_RB("LV",x,theta, obs_list_noise, lower_list, upper_list, tout_list, 1000, 50,diag(3))
covar2


