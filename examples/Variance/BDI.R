source("examples/Parameters/BDI/Equilibrium (N=1000).R")

tau = 0.1
lower_list = list()
upper_list = list()
inc = c(1)
n_inc = 61
obs = 1000
for(i in 0:(n_inc-1)){
  lower_list = append(lower_list,obs-inc*i)
  upper_list = append(upper_list,obs+inc*i)
}
taus = seq(0,1,by=0.01)

pdf("Test.pdf")
contour_data = get_contour_brownian("BDI", theta, tout, taus, x, obs, lower_list, upper_list)
filled.contour(contour_data,
plot.axes = {
  axis(1)
  axis(2)
  contour(contour_data,add=T,lwd=2,levels=c(seq(0,0.9,by=0.1),seq(0.91,1,by=0.01)))
})

contour_data = get_contour_brownian_fast("BDI", theta, tout, taus, x, obs, lower_list, upper_list)
filled.contour(contour_data,
plot.axes = {
  axis(1)
  axis(2)
  contour(contour_data,add=T,lwd=2,levels=c(seq(0,0.9,by=0.1),seq(0.91,1,by=0.01)))
})
dev.off()

microbenchmark(
  {get_contour_brownian("BDI", theta, tout, taus, x, obs, lower_list, upper_list)},
  {get_contour_brownian_fast("BDI", theta, tout, taus, x, obs, lower_list, upper_list)}
)

lower
upper
obs = 1005


get_variance_brownian("BDI",theta,tout,tau,x,obs,lower,upper)
get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,lower,upper)
tau = 0.01


