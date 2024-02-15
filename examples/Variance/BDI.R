source("examples/Parameters/BDI/Equilibrium (N=1000).R")

tau = 0.1
lower_list = list()
upper_list = list()
inc = c(1)
n_inc = 61
for(i in 0:(n_inc-1)){
  lower_list = append(lower_list,obs-inc*i)
  upper_list = append(upper_list,obs+inc*i)
}
taus = seq(0,1,by=0.01)

pdf("variance_contours.pdf")
lower_list
contour_data = get_contour_exact("BDI", theta, tout, taus, x, obs, lower_list, upper_list,lower_limit,upper_limit)
filled.contour(contour_data,
plot.axes = {
  axis(1)
  axis(2)
  contour(contour_data,add=T,lwd=2,levels=c(seq(0,0.9,by=0.1),seq(0.91,1,by=0.01)))
  title(main="Exact",cex.main=1.5,xlab="tau",ylab="widths",cex.lab=1.5)
})

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
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"PVR"=1-var_fast,"goal"=lvls[i],"tau"=lvls[k],"col"=col[k],"type"="brownian_fast"))
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"PVR"=1-var,"goal"=lvls[i],"tau"=lvls[k],"col"=col[k],"type"="brownian"))
    box_results = rbind(box_results,data.frame("lower"=results$lower,"upper"=results$upper,"PVR"=1-var_exact,"goal"=lvls[i],"tau"=lvls[k],"col"=col[k],"type"="exact"))
  }
  # lines(box_results$goal,box_results$PVR,col=col[k])
}

subset(box_results,goal==0.1)
box = subset(box_results,(tau==0.1 & type=="brownian"))
lines(box$goal,box$PVR,type='o',col="red",lwd=2)
box = subset(box_results,(tau==0.1 & type=="brownian_fast"))
lines(box$goal,box$PVR,type='o',col="blue",lwd=2)
box = subset(box_results,(tau==0.1 & type=="exact"))
lines(box$goal,box$PVR,type='o',col="green",lwd=2)
box
title(main=paste("N:",1000),cex.main=1.5)
legend("topleft",legend=c("BM sum over x_{t-tau}in S","BM Integrated out x_{t-tau}in Omega","Exact"),col=c("red","blue","green"),inset=0.05,lty=c(1,1),cex=1.5)
dev.off()

subset(box_results,tau==0.1)

get_variance_brownian_fast("BDI",theta,tout,0.1,x,obs,lower,upper)
get_variance_brownian("BDI",theta,tout,0.1,x,obs,lower,upper)
lower_limit = 960
upper_limit = 1040

get_variance_exact("BDI",theta,tout,0.1,x,obs,lower,upper,lower_limit,upper_limit)
upper

get_mu("BDI",x,theta)
