library(mvtnorm)
library(coda)

rwm=function(nits,theta0,log.pi,lambda,qV,...) {
  d=length(theta0)
  if (is.null(qV)) {
    A=diag(d)
  }
  else {
    ## If d=1, for Z~N(0,1), X=mu+aZ ~ N(mu,a^2).
    ## If d>1, for Z~MVN(0,I_d), X=mu+AZ ~ N(mu,AA').
    ## The Cholesky decomposition finds A such that A'A=qV.
    A=t(chol(qV))
  }
  store=matrix(nrow=nits+1,ncol=d+1)
  psis=matrix(nrow=nits,ncol=d)
  nacc=0
  theta.curr=theta0
  log.pi.curr=log.pi(theta.curr,...)
  store[1,]=c(theta.curr,log.pi.curr)
  for (i in 1:nits) {
    psi=theta.curr+lambda*A%*%rnorm(d); psis[i,]=psi
    log.pi.prop=log.pi(psi,...)
    log.alpha=log.pi.prop-log.pi.curr
    if (log(runif(1))<log.alpha) {
      theta.curr=psi
      log.pi.curr=log.pi.prop
      nacc=nacc+1
    }
    store[i+1,]=c(theta.curr,log.pi.curr)
  }
  return(list(acc=nacc/nits,store=store,psis=psis))
}

get_h = function(x,thetas){
  return(c(thetas[1]*x[1],thetas[2]*x[1]*x[2],thetas[3]*x[2]))
}
# gillespie_alg_entire = function(x0,time_step,thetas,S){
#   h = get_h(x0,thetas)  
#   h0 = sum(h)
#   probs = h/h0
#   curr_time = 0
#   next_reaction_time = curr_time + rexp(1,rate=h0)
#   counter = 1
#   xt = data.frame(matrix(nrow=0,ncol=nrow(S)+1))
#   xt[counter,] = c(0,x0) 
#   while(next_reaction_time<time_step){
#     curr_time = next_reaction_time
#     reaction = sample(1:3,1,replace=T,probs)
#     counter = counter + 1
#     xt[counter,1] = curr_time
#     xt[counter,-1] = xt[counter-1,-1] + S[,reaction]
#     next_reaction_time = curr_time + rexp(1,rate=h0)
#     h = get_h(as.numeric(xt[counter,-1]),thetas)  
#     h0 = sum(h)
#     probs = h/h0
#   }
#   return(xt)
# }

gillespie_alg_ = function(x0,time_step,thetas,S){
  if(sum(abs(x0-c(0,0)))==0){
    print(paste(xt,thetas,counter))
    return(x0)
  }
  h = get_h(x0,thetas)  
  h0 = sum(h)
  probs = h/h0
  curr_time = 0
  next_reaction_time = curr_time + rexp(1,rate=h0)
  counter = 1
  xt = x0
  while(next_reaction_time<time_step){
    if(sum(abs(xt-c(0,0)))==0 || counter>1000){
      return(xt)
    }
    curr_time = next_reaction_time
    reaction = sample(1:3,1,replace=T,probs)
    counter = counter + 1
    xt = xt + S[,reaction]
    next_reaction_time = curr_time + rexp(1,rate=h0)
    h = get_h(xt,thetas)  
    h0 = sum(h)
    probs = h/h0
  }
  return(xt)
}
log.pi(log(thetas))
log.pi((thetas))

thetas= c(0.5,0.0025,0.3)
S = matrix(c(1,-1,0,
             0,1,-1),nrow=2,ncol=3,byrow=T)
S = matrix(c(1,-1,1),ncol=3)
simX0 = function(thetas) {
  return(c(100,100))
}

simNextX = function(x,thetas) {
  xnew = gillespie_alg_(x,1,thetas,S)
  return(xnew)
}
simY = function(x,thetas) {
  num_species = nrow(S)
  return(x + rnorm(n = num_species,sd=NOISE))
}
logLik = function(y,x,thetas) {
  return(sum(dnorm(y,x,log=T,sd=NOISE)))
}

get_y = function(thetas,num_obs){
  num_species = nrow(S)
  x=simX0()
  y = matrix(ncol=num_species,nrow=num_obs)
  y[1,] = x+rnorm(num_species,sd=NOISE)
  for(i in 2:num_obs){
    x = simNextX(x,thetas)
    y[i,] = x+rnorm(num_species,sd=NOISE)
  }
  return(y)
}

BSfilter<-function(M,ys,simX0,simXt,loglik,...) {
  d=length(simX0(...)) ## dimension of state space
  logwstars=rep(0,M)
  nobs=nrow(ys)
  xs=matrix(nrow=M,ncol=d,data=0)
  lml=0 ## log marginal likelihood
  for (t in 1:nobs) {
    if(t==1) {
      for (i in 1:M) {
        xs[i,]=simX0(...)
      }
    }
    else {
      for (i in 1:M) {
        xs[i,]=simXt(xs[i,],...)
      }
    }
    for (i in 1:M) {
      logwstars[i]=loglik(ys[t,],xs[i,],...)
    }
    mws=max(logwstars)
    logwstars=logwstars-mws ## for robustness
    wstars=exp(logwstars)
    if(t>0){
      lml=lml+mws+log(mean(wstars))
    }
    ws=wstars/sum(wstars) ## normalised weights
    if(prod(is.nan(ws))==1){
      return(list(lml=-Inf,xs=NaN,ws=NaN))
      print("-Inf")
    }
    if (t<nobs) {
      ## get indices of resampled particles
      iresamp=sample(1:M,M,prob=ws,replace=TRUE)
      xsresamp=as.matrix(xs[iresamp,])
      xs=xsresamp
    }
  }
  return(list(lml=lml,xs=xs,ws=ws))
}


simNextX(simX0(),theta)
BSfilter(10,y,simX0,simNextX,logLik,thetas)$lml


sim_data = function(M,num_obs){
  data=list()
  data[[1]] = t(replicate(simX0(thetas),n=M))
  for(i in 2:num_obs){
    data[[i]] = t(vsimNextX(x1=data[[i-1]][,1],x2=data[[i-1]][,2],thetas))
  }
  return(data)
}
sim_y = function(data){
  num_obs = length(data)
  M = dim(data[[1]])[1]
  y = list()
  for(i in 1:num_obs){
    y_temp = matrix(0,nrow=M,ncol=1)
    for(m in 1:M){
      y_temp[m,] = simY(data[[i]][m,],thetas)
    }
    y[[i]] = y_temp
  }
  return(y)
}

log.pi = function(thetas){
  ll = BSfilter(M,y,simX0,simNextX,logLik,exp(thetas))$lml+
    sum(1/prod(exp(thetas))*dnorm(exp(thetas),mean = (c(0.5,0.0025,0.3)),sd = rep(1,3),log = T))
    return(ll)
  # 1/(|g'(\theta)|)f_X(g(\theta))
  # g(\theta)=\exp(\theta)
}
log.pi = function(thetas){
  ll = BSfilter(M,y,simX0,simNextX,logLik,exp(thetas))$lml+
    sum(dnorm(exp(thetas),mean = (c(0.5,0.0025,0.3)),sd = rep(2,3),log = T))
  return(ll)
  # 1/(|g'(\theta)|)f_X(g(\theta))
  # g(\theta)=\exp(\theta)
}
log.pi(log(thetas))

vsimNextX = Vectorize(function(x1,x2,thetas){simNextX(c(x1,x2),thetas)},vectorize.args = c('x1','x2'))
num_obs = 5
M = 500
#data = sim_data(M,num_obs)
NOISE = 2
y=get_y(thetas,num_obs)[-1,]
y
get_obs("LV",x,thetas,c(1,2,3,4))

z = rwm(nits = 50000,
        theta0 = log(thetas),
        log.pi = log.pi,
        lambda = 1,
        qV = qV
)
y = matrix(c(253,197,
             286,194,
             329,134,
             336,102))
y = matrix(c(254.8337,197.5254,
             284.6064,195.0731,
             332.5094,135.0130,
             339.3588,102.0116),nrow=4,ncol=2,byrow=T)
boot = (replicate(500,BSfilter(500,y,simX0,simNextX,logLik,thetas)$lml))
BSfilter(500,y,simX0,simNextX,logLik,thetas)$lml
y = obs_list + rnorm(4,sd=2)

var(boot)
mean(boot)
hist(boot)
num_obs=4

z$acc
plot(z$store[,3],type='l')
qV = cov(z$store[,1:3])

mean(exp(z$store[,1]))
mean(exp(z$store[,2]))
mean(exp(z$store[,3]))

quantile(exp(z$store[,1]),c(0.025,0.975))
quantile(exp(z$store[,2]),c(0.025,0.975))
quantile(exp(z$store[,1]),c(0.025,0.975))
thetas
thetas
effectiveSize(z$store[,1:3])