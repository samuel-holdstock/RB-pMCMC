source("examples/Parameters/BDI/Equilibrium (N=100).R")
tau = get_tau("BDI",x,theta,obs,tout,1000)
get_box_brownian_fast("BDI",theta,tout,tau,x,obs,2)

get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,86,114)
get_VRF_big_box_tau("BDI",x,theta,obs,tout,tau)
get_VRF_big_box_tau2("BDI",x,theta,obs,tout,tau)

get_tau_list("BDI",x,theta,obs_list,tout_list,sqrt(0.90))
get_box("BDI",theta,tout_list,tau,x,obs_list,sqrt(0.99))
get_PVR_big_box_tau("BDI",x,theta,obs,tout,1)
get_mu("BDI",x,theta)
get_var("BDI",x,theta)
get_covar("BDI",x,theta)


get_Q1_brownian_fast(1.28,theta,tout,tau,x,obs,284,36)
get_Q2_brownian_fast(1.28,theta,tout,tau,x,obs,110,284,36)
get_Q3_brownian_fast(1.28,theta,tout,tau,x,obs,100,284,36)

tau = 0.1
lower_list = list()
upper_list = list()
inc = c(1)
n_inc = 13
for(i in 0:(n_inc-1)){
  lower_list = append(lower_list,obs-inc*i)
  upper_list = append(upper_list,obs+inc*i)
}
taus = seq(0,1,by=0.001)
# Bayesian bootstrap
# Steve walker
# Generalised Bayes 

pdf("variance_contours.pdf")
contour_data = get_contour_exact("BDI", theta, tout, taus, x, obs, lower_list, upper_list,lower_limit,upper_limit)
contour_data = get_contour_brownian_fast("BDI", theta, tout, taus, x, obs, lower_list, upper_list)
X = dim(contour_data)[2]-1
Y = dim(contour_data)[1]
goals = seq(0,100,length.out=21)[-1]
goal = 10
alpha = seq(goal,goal+200,length.out=200)

png("variance_contours.png",res=600,width=10,height=10,units='in')
pdf("VRF_contours.pdf")
filled.contour(y=2*(0:X),x=seq(0,1,length.out=Y),z=contour_data,
plot.axes = {
  axis(1,cex.axis=2)
  axis(2,cex.axis=2)
  contour(y=2*(0:X),x=seq(0,1,length.out=Y),z=contour_data,add=T,lwd=2,levels=goals,labcex=2)
  title(main="Variance reduction factor (VRF)",cex.main=2,xlab="tau",ylab="widths",cex.lab=2,cex.main=2)
  L = 2000
  max_goal = goal+2000
  x_pos = rep(0,L)
  y_pos = rep(0,L)
  counter = 1
  for(alpha in seq(goal,max_goal,length.out=L)){
    pos=get_width_tau_point(goal,alpha)
    x_pos[counter] = pos[1]
    y_pos[counter] = pos[2]
    counter = counter+1
    # points(x=pos[1],y=pos[2],cex=2,col='cyan',pch=3,lwd=2)
  }
  # points(x=x_pos,y=y_pos,cex=2,col=c(rep('red',500),rep('blue',500)),lwd=2)
  # points(x=x_pos,y=y_pos,cex=2,col=color.gradient(seq(goal,max_goal,length.out=L)),lwd=2)
  segments(x0=head(x_pos,-1),x1=tail(x_pos,-1),y0=head(y_pos,-1),y1=tail(y_pos,-1),cex=2,col=color.gradient(seq(goal,max_goal,length.out=L-1)),lwd=2)
  # lines(x=x_pos,y=y_pos,cex=2,col='cyan',lwd=2,type='l')
})
dev.off()

color.gradient <- function(x, colors=c("green", "blue"), 
                           colsteps=100) {
  return(colorRampPalette(colors)(colsteps)[
    findInterval(x, seq(min(x), max(x), length.out=colsteps))
  ])
}

goal = 10
# Fix PVR
# Choose alpha
# That gives a beta
# This gives a tau
# This gives a width
# Plot (tau,width)

get_width_tau_point = function(goal,alpha){
  a = alpha
  b = alpha/(goal)
  tau = get_tau("BDI", x, theta, obs, tout, a)
  width = 2*(obs-get_box_brownian_fast("BDI",theta,tout,tau,x,obs,b)$lower)
  return(c(tau,width))
}
get_width_tau_point(goal,1000)

dev.off()
get_line = function(goal){
  tau = get_tau("BDI", x, theta, obs, tout, sqrt(goal))
  width = 2*(100-get_box_brownian_fast("BDI", theta, tout, tau, x, obs, sqrt(goal))$lower)
  return(c(tau,width))
}
get_line = function(alpha,goal){
  tau = get_tau("BDI", x, theta, obs, tout, alpha)
  beta = alpha/goal
  width = 2*(100-get_box_brownian_fast("BDI", theta, tout, tau, x, obs, beta)$lower)
  return(c(tau,width))
}
tau=get_tau("BDI",x,theta,obs,tout,100)
get_tau("BDI",x,theta,obs,tout,10000)

get_box_brownian_fast("BDI", theta, tout, tau, x, obs,5)
get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,96,104)
get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,95,105)

vget_variance_brownian_fast = Vectorize(get_variance_brownian_fast,vectorize.args = c('tau'))
vget_Q2_brownian_fast = Vectorize(get_Q2_brownian_fast,vectorize.args = c('upper'))
vget_Q3_brownian_fast = Vectorize(get_Q3_brownian_fast,vectorize.args = c('lower'))

vget_variance_brownian_fast("BDI",theta,tout,0.8459,x,obs,30,180)
RB("BDI",x,theta,tout,30,180,0.94,obs,10)

p = RB("BDI",x,theta,tout,30,180,1,obs,10)
plot(taus,vget_variance_brownian_fast("BDI",theta,tout,taus,x,obs,50,150),type='l',ylab='Variance of infintely wide box (normalised)',xlab='Tau',cex.lab=1.5,cex.axis=1.5)
plot(seq(0,0.95,length.out=100),vget_variance_brownian_fast("BDI",theta,tout,seq(0,0.8,length.out=100),x,obs,50,150),type='l',ylab='Variance of infintely wide box (normalised)',xlab='Tau',cex.lab=1.5,cex.axis=1.5)
abline(h=0.05,lty='dashed')
mu = get_mu("BDI",x,theta)
sig2 = get_var("BDI",x,theta)
widths = seq(2,20,by=2)
Q1 = get_Q1_brownian_fast(1.1,theta,tout,tau,x,obs,mu,sig2)
Q2 = vget_Q2_brownian_fast(1.1,theta,tout,tau,x,obs,x+(0:10),mu,sig2)
Q3 = vget_Q3_brownian_fast(1.1,theta,tout,tau,x,obs,x-(0:10),mu,sig2)
plot(2*(0:10),(Q2+Q3),type='o',xlab='width',ylab='Diff variance',cex.lab=1.5,cex.axis=1.5)
abline(h=(beta-1)/beta*(p*(1-p)-Q1),lty='dashed')
beta = 10
# xaxis alpha, yaxis tau

plot_alpha_tau = function(goal,alpha){
  a = alpha
  b = a/goal
  # b = seq(0,1/(1-goal),length.out=999)[-1]
  # b = b[(1-goal)*b<1-goal]
  # a = (1-goal)*b
  vget_tau = Vectorize(get_tau,vectorize.args = c('alpha'))
  tau = vget_tau("BDI",x,theta,obs,tout,a)
  plot(a,tau,type='l',xlab='alpha',ylab='tau')
  # plot(b,tau,type='l',xlab='beta',ylab='tau')
}
plot_alpha_tau(goal,alpha)

# xaxis alpha, yaxis width
goal = 0.9
plot_alpha_width = function(goal,alpha){
  a = alpha
  b = a/goal
  # b = seq(0,1,length.out=999)[-c(1,999)]
  # b = seq(0,1/(1-goal),length.out=999)[-1]
  # b = b[(1-goal)*b<1-goal]
  # a = (1-goal)*b
  vget_tau = Vectorize(get_tau,vectorize.args = c('alpha'))
  tau = vget_tau("BDI",x,theta,obs,tout,a)
  widths = 2*(obs-as.numeric(mapply(get_box_brownian_fast,"BDI",list(theta),tout,tau,x,obs,b)[1,]))
  plot(a,widths,type='o',xlab='alpha',ylab='width')
  # plot(b,widths,type='o',xlab='beta',ylab='width')
}
plot_alpha_width(goal,alpha)
tau=get_tau("BDI",x,theta,obs,tout,10)
get_box_brownian_fast("BDI",theta,tout,tau,x,obs,1)

# xaxis alpha, yaxis width
goal = 0.95
plot_alpha_PVR = function(goal){
  b = seq(0,1,length.out=999)[-c(1,999)]
  b = b[(1-goal)*b<1-goal]
  a = (1-goal)*b
  vget_tau = Vectorize(get_tau,vectorize.args = c('goal'))
  tau = vget_tau("BDI",x,theta,obs,tout,1-a)
  widths = 2*(obs-as.numeric(mapply(get_box_brownian_fast,"BDI",list(theta),tout,tau,x,obs,b)[1,]))
  lower = obs-widths/2
  upper = obs+widths/2
  PVR = mapply(get_variance_brownian_fast,"BDI",list(theta),tout,tau,x,obs,lower,upper)
  plot(a,1-PVR,type='o',xlab='alpha',ylab='width')
  # plot(b,1-PVR,type='o',xlab='beta',ylab='width')
}
plot_alpha_PVR(0.99)
abline(h=goal)
get_variance_brownian_fast("BDI",theta,tout,0.1,x,obs,90,110)
obs=110

alpha = 0.1
beta = 0.95
1-alpha/beta

(1-alpha) # reduce by
(1-beta) # increase by



goals = seq(0.1,0.9,length.out=9)
goals = seq(0.9,1,length.out=10)[-10]
goals = c(seq(0.1,0.9,length.out=9),seq(0.9,1,length.out=10)[-10])
for(i in 1:length(goals)){
  pos = get_line(goals[i])
  points(x=pos[1],y=pos[2],cex=3,col='cyan',pch=4,lwd=3)
}


contour_data = get_contour_brownian("BDI", theta, tout, taus, x, obs, lower_list, upper_list)
filled.contour(contour_data,
plot.axes = {
  axis(1)
  axis(2)
  contour(contour_data,add=T,lwd=2,levels=c(seq(0,0.9,by=0.1),seq(0.91,1,by=0.01)))
  title(main="Brownian sum",cex.main=1.5,xlab="tau",ylab="widths",cex.lab=1.5)
})

contour_data = get_contour_brownian_fast("BDI", theta, tout, taus, x, obs, lower_list, upper_list)
filled.contour(contour_data,
plot.axes = {
  axis(1)
  axis(2)
  contour(contour_data,add=T,lwd=2,levels=c(seq(0,0.9,by=0.1),seq(0.91,1,by=0.01)))
  title(main="Brownian integral",cex.main=1.5,xlab="tau",ylab="widths",cex.lab=1.5)
})
dev.off()


# microbenchmark(
#   {get_contour_brownian("BDI", theta, tout, taus, x, obs, lower_list, upper_list)},
#   {get_contour_brownian_fast("BDI", theta, tout, taus, x, obs, lower_list, upper_list)},
#   {get_contour_exact("BDI", theta, tout, taus, x, obs, lower_list, upper_list,lower_limit,upper_limit)},unit='micro'
# )

lower = c(990)
upper = c(1010)

get_variance_brownian("BDI",theta,tout,tau,x,obs,lower,upper)
get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,lower,upper)

mu = get_mu("BDI",x,theta)
var = get_var("BDI",x,theta)
get_Q1_brownian_fast("BDI",theta,tout,tau,x,obs,mu,var)
get_Q2_brownian_fast("BDI",theta,tout,tau,x,obs,obs+5,mu,var)
get_Q3_brownian_fast("BDI",theta,tout,tau,x,obs,obs-5,mu,var)

box_results = list()
lvls=seq(0,0.9,by=0.1)
col=1:length(lvls)
for(tau in seq(0.1,1,by=0.1)){
  for(i in 1:length(lvls)){
    results = get_box_brownian_fast("BDI",theta,tout,tau,x,obs,lvls[i])
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"tau"=tau,"col"=col[i]))
  }
}
xvals = split(box_results$tau,box_results$col)
yvals = split((box_results$upper-obs+1)/61,box_results$col)
contour_data = get_contour_brownian_fast("BDI", theta, tout, taus, x, obs, lower_list, upper_list)
filled.contour(contour_data,
plot.axes = {
  axis(1)
  axis(2)
  contour(contour_data,add=T,lwd=2,levels=c(seq(0,0.9,by=0.1),seq(0.91,1,by=0.01)))
  points(box_results$tau,(box_results$upper-obs+1)/61,col=box_results$col,pch=4,cex=2,lwd=2)
  mapply(lines, xvals, yvals,col=box_results$col)
}
)

source("examples/Parameters/BDI/Equilibrium (N=1000).R")
pdf("Test.pdf")

plot(x=-1,xlim=c(0,1),ylim=c(0,1),xlab="Percentage reduction in variance target",ylab="Percentage reduction in variance acheived",cex.lab=1.5,cex.axis=1.5)
lines(x=c(0,1),y=c(0,1),lty='dashed')
box_results = data.frame()
for(k in 2:length(lvls)){
  for(i in 1:length(lvls)){
    results = get_box_brownian_fast("BDI",theta,tout,lvls[k],x,obs,lvls[i])
    var_fast = get_variance_brownian_fast("BDI",theta,tout,lvls[k],x,obs,results$lower,results$upper)
    var = get_variance_brownian("BDI",theta,tout,lvls[k],x,obs,results$lower,results$upper)
    var_exact = get_variance_exact("BDI",theta,tout,lvls[k],x,obs,results$lower,results$upper,lower_limit,upper_limit) 
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"PVR"=1-var_fast,"goal"=lvls[i],"scale"=results$scale,"tau"=lvls[k],"col"=col[k],"type"="brownian_fast"))
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"PVR"=1-var,"goal"=lvls[i],"scale"=results$scale,"tau"=lvls[k],"col"=col[k],"type"="brownian"))
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"PVR"=1-var_exact,"goal"=lvls[i],"scale"=results$scale,"tau"=lvls[k],"col"=col[k],"type"="exact"))
  }
  # lines(box_results$goal,box_results$PVR,col=col[k])
}

box = subset(box_results,(tau==0.1 & type=="brownian"))
lines(box$goal,box$PVR,type='o',col="red",lwd=2)
box = subset(box_results,(tau==0.1 & type=="brownian_fast"))
lines(box$goal,box$PVR,type='o',col="blue",lwd=2)
box = subset(box_results,(tau==0.1 & type=="exact"))
lines(box$goal,box$PVR,type='o',col="green",lwd=2)

lines(box$goal,box$PVR*box$scale,type='o',col="red",lwd=2)
box = subset(box_results,(tau==0.1 & type=="brownian_fast"))
lines(box$goal,box$PVR*box$scale,type='o',col="blue",lwd=2)
box = subset(box_results,(tau==0.1 & type=="exact"))
lines(box$goal,box$PVR*box$scale,type='o',col="green",lwd=2)

get_mu("BDI",x,theta)
get_var("BDI",x,theta)
get_rate("BDI",x,theta)
1/2+1/2*0.8*sqrt(get_var("BDI",x,theta))/sqrt(sum(get_rate("BDI",x,theta)))
1/2+1/2*0.8*sqrt(400)/sqrt(400)

box
title(main=paste("N:",1000),cex.main=1.5)
legend("topleft",legend=c("BM sum over x_{t-tau}in S","BM Integrated out x_{t-tau}in Omega","Exact"),col=c("red","blue","green"),inset=0.05,lty=c(1,1),cex=1.5)
dev.off()

subset(box_results,tau==0.1)

get_variance_brownian_fast("SIR_N",theta,tout,tau,x,obs,lower_limit,upper_limit)

get_variance_brownian("SIR_N",theta,tout,0.1,x,obs,lower,upper)
lower_limit = 960
upper_limit = 1040

get_variance_exact("SIR_N",theta,tout,0.1,x,obs,lower,upper,lower_limit,upper_limit)
get_box_brownian_fast("SIR_N",theta,tout,tau,x,obs,0.1)

source("examples/Parameters/SIR_N/Equilibrium (N=1000).R")

qnorm(1-0.025/2)

