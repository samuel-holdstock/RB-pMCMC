theta = c(0.3,0.2)
mu = get_mu("BDI_2",x,theta)
sig = sqrt(get_var("BDI_2",x,theta))

lower = 950
upper = 1010
s = seq(lower,upper,length.out=100)
r = seq(-5,40,length.out=100)
lower_limit = 920
upper_limit = 1080

((1+mu^2/sig^4)*pnorm(mu/sig^2)+mu/(sig^2*sqrt(2*pi))*exp(-1/2*mu^2/(sig^4)))/(2/sqrt(2*pi)*exp(-mu^2/(2*sig^4))+2*mu/sig^2*pnorm(mu/sig^2))+0.5
sqrt(2*pi)/4+0.5
(1/2+2*mu/sqrt(2*pi*sig^4))/(2/sqrt(2*pi)+2*mu/sig^2)+0.5
mu_h = mu/sig^2
((1+mu_h^2)*pnorm(mu_h)+mu_h/sqrt(2*pi)*exp(-1/2*mu_h^2))/(2*(1/sqrt(2*pi)*exp(-1/2*mu_h^2)+mu_h*pnorm(mu_h)))+0.5


get_estimate_continuous = function(m,wt,start,mu,sig,tau,offset){
  vget_above_maximum_sigma_wt(m-start+offset,mu,sig,wt-start,tau)
}
get_estimate_discrete = function(m,wt,start,mu,sig,tau){
  if(wt>m || start>m){
    return(expm(get_coffin_matrix("BDI_2",lower_limit,upper_limit,theta)*tau)[start-lower_limit+1,wt-lower_limit+1])
  }
  (expm(get_coffin_matrix("BDI_2",lower_limit,upper_limit,theta)*tau)[start-lower_limit+1,wt-lower_limit+1])-
  (expm(get_coffin_matrix("BDI_2",lower_limit,m,theta)*tau)[start-lower_limit+1,wt-lower_limit+1])
}
1001-lower_limit+1

get_estimate_continuous2 = function(m,wt,start,mu,sig,tau){
  p = ((mu/sig^2)+1)/2
  q = 1-p
  n = sig^2*tau
  m = m-start
  wt = wt-start
  (p/q)^m * 1/sqrt(n*p*q)*dnorm((n+wt-2*m-2*n*p)/(2*sqrt(n*p*q)))/2
}
get_estimate_continuous(upper,obs,x,mu,sig,tau,1) # P(M>=u+1, X_t=obs| X_0=x)
get_estimate_continuous2(upper,obs,x,mu,sig,tau) # P(M>=u, X_t=obs| X_0=x)
get_estimate_discrete(upper,obs,x,mu,sig,tau) # P(M>=u+1, X_t=obs| X_0=x) = P(M<=u, X_t=obs| X_0=x) - P(M>u, X_t=obs| X_0=x)

get_table = function(ms){
  results = data.frame(matrix(ncol=3,nrow=0))
  colnames(results) = c('disc','cont','cont_offset')
  for(i in (1:length(ms))){
    m = ms[i]
    results[i,] = c(get_estimate_discrete(m,wt,x,mu,sig,tau),
    get_estimate_continuous(m,wt,x,mu,sig,tau,0),
    get_estimate_continuous(m,wt,x,mu,sig,tau,0))
  }
  results
}
results = get_table(ms)
plot(results$cont,results$disc,type='o',col='blue')
lines(results$cont_offset,results$disc,type='o',col='red')
abline(a=0,b=1,lty='dashed')

ms
tau = 2
tau=0.1
get_estimate_continuous(1050,wt,x,mu,sig,tau,0.13)
get_estimate_continuous(1050,wt,x,mu,sig,tau,-0.13)

vget_estimate_discrete = Vectorize(get_estimate_discrete,vectorize.args = 'm')
vget_estimate_continuous = Vectorize(get_estimate_continuous,vectorize.args = 'm')
vget_estimate_continuous2 = Vectorize(get_estimate_continuous2,vectorize.args = 'm')

tau=0.5
obs = round(x + mu*tau)
ms = obs+0:30

plot(ms,log(vget_estimate_discrete(ms,obs,x,mu,sig,tau),10),type='o')
((1+mu_h^2)*pnorm(mu_h)+mu_h/sqrt(2*pi)*exp(-1/2*mu_h^2))/(2*(1/sqrt(2*pi)*exp(-1/2*mu_h^2)+mu_h*pnorm(mu_h)))+0.5
lines(ms,log(vget_estimate_continuous(ms,obs,x,mu,sig,tau,1.171381),10),col='red')
lines(ms,log(vget_estimate_continuous(ms,obs,x,mu,sig,tau,0.5),10),col='green')
box = get_box_brownian_fast("BDI_2",theta,tout,tau,x,obs,0.90)
abline(v=box$upper)

lines(ms,log(vget_estimate_continuous(ms,obs,x,mu,sig,tau,0.5),10),col='cyan')
lines(ms,vget_estimate_continuous(ms,obs,x,mu,sig,tau,0),col='blue')
lines(ms,vget_estimate_continuous2(ms+0.57,obs,x,mu,sig,tau),col='darkgreen')
lines(ms,vget_estimate_continuous2(ms+1,obs,x,mu,sig,tau),col='orange')
lines(ms,vget_estimate_continuous2(ms,obs,x,mu,sig,tau),col='purple')

tau=1

upper_limit = 1260
lower_limit = 950
plot(ms,vget_estimate_discrete(ms,obs,x,mu,sig,tau),type='o')
lines(ms,vget_estimate_continuous(ms,obs,x,mu,sig,tau,1),col='purple')

tau = 2
get_above_maximum_wt2 = function(x,mu,sig,wt,tout){
  if(tout<=1e-15){
    return(0)
  }
  k = max(wt,0,x)
  probability = sqrt(1/(2*pi*sig^2*tout))*exp(mu*wt/sig^2 - mu*mu*tout/(2*sig^2) - (2*k-wt)*(2*k-wt)/(2*tout*sig^2))
  return(probability)
}
get_above_maximum_wt3 = function(m,mu,sig,wt,tout){
  if(tout<=1e-15){
    return(0)
  }
  probability = exp(2*mu*m/(sig^2))*dnorm((wt-2*m-mu*tout)/(sig*sqrt(tout)))/(sig*sqrt(tout))
  return(probability)
}
get_above_maximum_wt2(1045-start,mu,sig,wt-start,tau)
get_above_maximum_wt3(1045-start,mu,sig,wt-start,tau)
tau = 0.1

wt = 1040
obs = 1040

start
wt
