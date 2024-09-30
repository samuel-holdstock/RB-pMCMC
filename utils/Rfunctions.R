get_chain = function(data,endTime){
    times = data[,1]    
    times = rep(times,each=2)[-1]
    times[length(times)+1] = endTime

    result = matrix(nrow=nrow(data)*2,ncol=ncol(data))
    for(i in 2:ncol(data)){
        result[,i] = rep(data[,i],each=2)
    }
    result[,1] = times
    return(result)
}

plot_chain = function(chain,species,colours,...){
    plot(chain[,c(1,2)],col=colours[1],type='l',...)
    columns = ncol(chain)
    if(columns>2){
        for(i in 3:columns){
            lines(chain[,c(1,i)],col=colours[i-1])
        }
    }
    legend("topleft",legend=species,col=colours,inset=0.05,lty=rep(1,columns-1))
}

plot_chain_list = function(chain_list,species,colours,...){
    plot(1,xlim=c(0,1),ylim=c(900,1100))
    for(chain in chain_list){
        columns = ncol(chain)
        for(i in 2:columns){
            lines(chain[,c(1,i)],col=colours[i-1])
        }
        legend("topleft",legend=species,col=colours,inset=0.05,lty=rep(1,columns-1))
    }
}

add_chain = function(chain,species,colours,...){
    lines(chain[,c(1,2)],col=colours[1],...)
    columns = ncol(chain)
    if(columns>2){
        for(i in 3:columns){
            lines(chain[,c(1,i)],col=colours[i-1])
        }
    }
}

add_box = function(lower,upper,tout,tau,colours,...){
    nboxes = length(lower)
    for(i in 1:nboxes){
        lines(c(tout,tout-tau,tout-tau,tout),c(upper[i],upper[i],lower[i],lower[i]),col=colours[i],lwd=3)
    }
}

add_box_goal = function(model_name,goal,theta,x,obs,tout,tau,colours){
    results = get_box_brownian_fast(model_name,theta,tout,tau,x,obs,goal)
    lower = results$lower
    upper = results$upper
    lines(c(tout,tout-tau,tout-tau,tout),c(upper,upper,lower,lower),col=colours,lwd=1.5)
}

add_obs = function(tout,obs,colours){
    npoints = length(obs)
    for(i in 1:npoints){
        points(tout,obs[i],col=colours[i],pch=4,cex=3)
    }
}

f = function(s, thetas, tout, tau, start, target, lower, upper){
    mu = get_mu("BDI",start,thetas)
    sig = sqrt(get_var("BDI",start,thetas))
    a = dnorm(s-start,mu*(tout-tau),sig*sqrt(tout-tau))
    b = dnorm(target-s,mu*tau,sig*sqrt(tau))   
    k = max(target, upper+1)
    v = min(target, lower-1)
    C = exp(-2*(k-target)*mu/(sig^2))
    alpha = dnorm(s-(target-mu*tau),0,sig*sqrt(tau))*(s>=k) + dnorm(s-(2*k-target-mu*tau),0,sig*sqrt(tau))*C*(s<k)
    D = exp(-2*(v-target)*mu/(sig^2))
    beta = dnorm(s-(target-mu*tau),0,sig*sqrt(tau))*(s<=v) + dnorm(s-(2*v-target-mu*tau),0,sig*sqrt(tau))*D*(s>v)
    delta = alpha + beta
    # a_ = upper-s+1
    # b_ = target-s
    # alpha = vget_above_maximum_sigma_wt(a_,mu,sig,b_,tau)
    # a_ = lower-s-1
    # b_ = target-s
    # beta = vget_below_minimum_sigma_wt(a_,mu,sig,b_,tau)
    # a = a*(a<1)+1*(a>=1)
    # b = b*(b<1)+1*(b>=1)
    # alpha = alpha*(alpha<1)+1*(alpha>=1)
    # beta = beta*(beta<1)+1*(beta>=1)
    # delta = delta*(delta<1)+1*(delta>=1)
    Q1 = a*b*(1-b)
    Q2 = a*alpha*(1-alpha)
    Q3 = a*beta*(1-beta)
    Q4 = a*delta*(1-delta)
    p = pnorm(target+1/2-start,mu*tout,sig*sqrt(tout))-pnorm(target-1/2-start,mu*tout,sig*sqrt(tout))
    percentage_variance_reduction =  (Q1-Q2-Q3)/(p*(1-p))
    # percentage_variance_reduction =  (Q1-Q4)/(p*(1-p))
    # return(Q2)
    return(percentage_variance_reduction)
}

vget_above_maximum_sigma_wt = Vectorize(get_above_maximum_sigma_wt,vectorize.args = c('x','wt'))
vget_below_minimum_sigma_wt = Vectorize(get_below_minimum_sigma_wt,vectorize.args = c('x','wt'))
vf = Vectorize(f,vectorize.args = 's')
obs = 1000
lower = 985
upper = 1055

f(1000,thetas=theta,tout=tout,tau=tau,start=x,target=obs,lower,upper)
f2(1000,thetas=theta,tout=tout,tau=tau,start=x,target=obs,lower,upper)

1-integrate(f,obs-100,obs+100,thetas=theta,tout=tout,tau=tau,start=x,target=obs,lower,upper,abs.tol = 1e-20)$value
# 1-integrate(f2,obs-100,obs+100,thetas=theta,tout=tout,tau=tau,start=x,target=obs,lower,upper,abs.tol = 1e-20)$value

get_variance_brownian_fast("BDI",theta,tout,tau,x,obs,lower,upper)
get_variance_brownian("BDI",theta,tout,tau,x,obs,lower,upper)
get_variance_exact("BDI",theta,tout,tau,x,obs,lower,upper,lower_limit,upper_limit)
lower = 995
upper = 1005

tau = 0.1

lower = 990
upper = 1020
obs = 1000
s = seq(lower-5,upper+5,by=0.01)

plot(s,vf(s,thetas=theta,tout=tout,tau=tau,start=x,target=obs,lower,upper),type='l')
lines(s,vf2(s,thetas=theta,tout=tout,tau=tau,start=x,target=obs,lower,upper),type='l',col='blue')
df = get_variance_brownian2("BDI",theta,tout,tau,x,obs,lower,upper)
df2 = get_variance_exact2("BDI",theta,tout,tau,x,obs,lower,upper,lower_limit,upper_limit)

points(lower:upper,(df$a*df$b*(1-df$b)-df$a*df$delta*(1-df$delta))/(df$p*(1-df$p)))
lines(type='o',lower:upper,(df2$a*df2$b*(1-df2$b)-df2$a*df2$delta*(1-df2$delta))/(df2$p*(1-df2$p)),col='darkgreen')
rect(xleft=lower:upper-0.5,ybottom=0,xright=lower:upper+0.5,ytop=(df2$a*df2$b*(1-df2$b)-df2$a*df2$delta*(1-df2$delta))/(df2$p*(1-df2$p)),col='darkgreen',density=0)
rect(xleft=lower:upper-0.5,ybottom=0,xright=lower:upper+0.5,ytop=(df$a*df$b*(1-df$b)-df$a*df$delta*(1-df$delta))/(df$p*(1-df$p)))
lower = 985
upper = 1015
vf(980,thetas=theta,tout=tout,tau=tau,start=x,target=obs,990,1020)
p = 0.01628562
df$pvr = (df$a*df$b*(1-df$b) - df$a*df$delta*(1-df$delta))/(p*(1-p))
plot(990:1020,df$pvr)

df$a*df$b*(1-df$b)
plot(990:1020,df$a*df$alpha*(1-df$alpha),type='o',col='blue')
lines(990:1020,vf(990:1020,thetas=theta,tout=tout,tau=tau,start=x,target=obs,990,1020),type='o',col='darkgreen')
legend("topleft",legend=c("Brownian density a_ialpha_i(1-alpha_i)","Brownian sum terms a_ialpha_i(1-alpha_i)"),col=c("blue","darkgreen"),inset=0.05,lty=c(1,1))

df$a*df$b*(1-df$b)
df$a*df$delta*(1-df$delta)
df$delta
vf(990:1020,thetas=theta,tout=tout,tau=tau,start=x,target=obs,990,1020)


n <- 400
t <- 1

sim_brown = function(m){
    bm <- x + c(0, cumsum(rnorm(n,mu*tout/n,sig*sqrt(tout/n))))
    # if(any(bm>m+1/2) & any(abs(bm[n+1] - wt)<1/2)){
    #     return(1)
    # }
    if(any(bm>=m)){
        return(1)
    }
    return(0)
}
sim_MJP = function(m,wt){
    upper = m
    df = sim_data("BDI_2",x,theta, tout, lower, upper, tau)
    if(!df$inS & any(tail(df$data,1)[,2]-wt==0)){
        return(1)
    }
    return(0)
}

B = 500000
sig = sqrt(get_var("BDI_2",x,theta))
mu = get_mu("BDI_2",x,theta)

sim1 = rep(NA,B)
for(i in 1:B){
    sim1[i] = sim_brown(m=upper)
}

sim2 = rep(NA,B)
for(i in 1:B){
    sim2[i] = sim_MJP(m=upper,wt=lower_limit:upper_limit)
}

sum(vexact_brown(upper-x+1,mu,sig,(lower_limit:upper_limit),tout))
sum(vexact_brown(upper-x+1/2,mu,sig,(lower_limit:upper_limit),tout))
sum(vexact_brown(upper-x,mu,sig,(lower_limit:upper_limit),tout))

sum(Vget_above_maximum_sigma_wt(upper-x+1,mu,sig,lower_limit:upper_limit,tout))
sum(Vget_above_maximum_sigma_wt(upper-x+1/2,mu,sig,lower_limit:upper_limit,tout))
sum(Vget_above_maximum_sigma_wt(upper-x,mu,sig,lower_limit:upper_limit,tout))
integrate(Vget_above_maximum_sigma_wt,lower_limit,upper_limit,x=upper-x+1,mu=mu,sig=sig,tout=tout,abs.tol=1e-15)$value

get_rate("BDI_2",x,theta)
theta = c(30,20)
upper = 140
mu = get_mu("BDI_2",x,theta)
sig = sqrt(get_var("BDI_2",x,theta))

mean(sim1) + c(-1,1)*1.96*sd(sim1)/sqrt(B)
mean(sim2) + c(-1,1)*1.96*sd(sim2)/sqrt(B)

max_test = function(upper,xT,x0,tout,mu,sig){
    1/sqrt(2*pi*sig^2*tout)*exp(-mu^2*tout/(2*sig^2))*exp(mu*(xT-x0)/sig^2)*exp(-((2*max(upper,x0,xT)-xT-x0)^2)/(2*sig^2*tout))
}
min_test = function(lower,xT,x0,tout,mu,sig){
    1/sqrt(2*pi*sig^2*tout)*exp(-mu^2*tout/(2*sig^2))*exp(mu*(xT-x0)/sig^2)*exp(-((2*min(lower,x0,xT)-xT-x0)^2)/(2*sig^2*tout))
}
sig = sqrt(get_var("BDI_2",x,theta))
mu = get_mu("BDI_2",x,theta)

max_test(3+0.5,3,0,1,mu,sig)
get_above_maximum_sigma_wt(3+0.5,mu,sig,3,1)

min_test(0,-3,0,1,mu,sig)
get_below_minimum_sigma_wt(0,mu,sig,-3,1)

vmax_text = Vectorize(max_test,vectorize.args = 'xT')
m=10
sum(max_test(1010+1/2,1000:1005,1000,tout,mu,sig))

lower = 900
upper = 1005
tau = 1

Q = get_coffin_matrix("BDI",lower_limit,upper_limit,theta)
tau=0.1
tail(c(expm(tau*Q)[,1]))
tail(c(vT_exp_Q(matrix(c(1,rep(0,nrow(Q)-1)),ncol=1),t(tau*Q),1e-40,F,T,F)))

get_Q = function(x,theta,obs,lower,upper,lower_limit,upper_limit){
    Q1 = get_coffin_matrix("BDI_2",lower_limit,upper,theta)
    Q2 = get_coffin_matrix("BDI_2",lower_limit,upper_limit,theta)
    bridge_index = state_to_index(upper,lower_limit,upper)+1
    target_index = state_to_index(obs,lower_limit,upper_limit)+1
    start_index = state_to_index(x,lower_limit,upper_limit)+1
    Q1[bridge_index,bridge_index+1]
    Q = adiag(Q1,Q2)
    Q1_n = nrow(Q1)
    Q2_n = nrow(Q2)
    Q_n = nrow(Q)
    bridge_index+1
    Q[bridge_index,bridge_index+1+bridge_index+1] = Q[bridge_index,bridge_index+1] 
    Q[bridge_index,bridge_index+1] = 0 
    return(list(expQ=expm(Q),Q=Q,bridge_index=bridge_index,target_index=target_index,start_index=start_index))
}
plot_Q = function(x,theta,obs,lower,upper,lower_limit,upper_limit){
    df = get_Q(x,theta,obs,lower,upper,lower_limit,upper_limit)
    expQ = df$expQ
    Q = df$Q
    bridge_index = df$bridge_index
    target_index = df$target_index
    Q_n = nrow(Q)
    mu = get_mu("BDI_2",x,theta)
    sig = sqrt(get_var("BDI_2",x,theta))   
    plot((upper+1):upper_limit,expQ[(bridge_index+1+bridge_index+1):(Q_n-1),bridge_index+1+target_index],type='l',xlim=c(lower_limit,upper_limit),col='blue',xlab="x_{t-tau}",ylab="Probability",sub="theta=c(220,200),t=tau=1,u=120, obs=110")
    # legend("topleft",c(
    #     "P(X_t=x_t,X_{[t-tau,t]}>u+1|X_{t-tau}=x_{t-tau}) for a BM",
    #     "P(X_t=x_t,X_{[t-tau,t]}>u|X_{t-tau}=x_{t-tau}) for a MJP"),col=c("blue","darkgreen"),lty=c(1,0),pch=c(NA,4))
    lines(lower_limit:upper,expQ[1:(bridge_index),bridge_index+1+target_index],type='l',col='blue')
    points(lower_limit:upper_limit,vget_above_maximum_sigma_wt(upper-(lower_limit:upper_limit)+1,mu,sig,obs-(lower_limit:upper_limit),1),col='darkgreen',pch=4)
}
upper = 120
obs = 110
pdf("Joint_brownian_max.pdf")
theta = c(220,200)
obs=95
plot_Q(x,theta,obs,lower,upper,lower_limit,upper_limit)
abline(v=upper)
dev.off()
lower_limit

get_Q = function(x,theta,obs,lower,upper,lower_limit,upper_limit){
    Q1 = get_coffin_matrix("BDI_2",lower_limit,upper,theta)
    Q2 = get_coffin_matrix("BDI_2",lower_limit,upper_limit,theta)
    target_index_Q2 = state_to_index(obs,lower_limit,upper_limit)+1
    start_index_Q2 = state_to_index(x,lower_limit,upper_limit)+1
    coffin_index_Q1 = upper-lower_limit+2
    coffin_index_Q2 = upper_limit-lower_limit+2
    qi = expm(Q1)
    bi = expm(Q2)
    qi = qi[-c(coffin_index_Q1),]
    bi = bi[-c(coffin_index_Q2),]
    qi = qi[,-c(coffin_index_Q1)]
    bi = bi[,-c(coffin_index_Q2)]
    qi = adiag(qi,diag(0,nrow(bi)-nrow(qi)))    
    di = bi-qi
    return(di)
}
plot_Q = function(x,theta,obs,lower,upper,lower_limit,upper_limit){
    expQ = get_Q(x,theta,obs,lower,upper,lower_limit,upper_limit)
    mu = get_mu("BDI_2",x,theta)
    sig = sqrt(get_var("BDI_2",x,theta))   
    obs_index = state_to_index(obs,lower_limit,upper_limit)+1
    plot(NULL,xlim=c(lower_limit,upper_limit),ylim=c(0,max(expQ[,obs_index])))
    legend("topleft",c(
        "P(X_t=x_t,X_{[t-tau,t]}>u+1|X_{t-tau}=x_{t-tau}) for a BM",
        "P(X_t=x_t,X_{[t-tau,t]}>u|X_{t-tau}=x_{t-tau}) for a MJP"),col=c("blue","darkgreen"),lty=c(1,0),pch=c(NA,4))
    lines(lower_limit:upper_limit,expQ[,obs_index],col='blue')
    points(lower_limit:upper_limit,vget_above_maximum_sigma_wt(upper-(lower_limit:upper_limit)+1,mu,sig,obs-(lower_limit:upper_limit),1),col='darkgreen',pch=4)
    points(lower_limit:upper_limit,vexact_brown(upper-(lower_limit:upper_limit)+1,mu,sig,obs-(lower_limit:upper_limit),1),col='darkred',pch=4)
}
theta = c(210,200)
plot_Q(x,theta,obs,lower,upper,lower_limit,upper_limit)
obs = 120
lower
upper = 121
x

get_Q_max = function(x,theta,obs,lower,upper,lower_limit,upper_limit){
    Q1 = get_coffin_matrix("BDI_2",lower_limit,upper,theta)
    Q2 = get_coffin_matrix("BDI_2",lower_limit,upper_limit,theta)
    start_index_Q1 = state_to_index(x,lower_limit,upper)+1
    start_index_Q2 = state_to_index(x,lower_limit,upper_limit)+1
    obs_index_Q1 = state_to_index(obs,lower_limit,upper)+1
    obs_index_Q2 = state_to_index(obs,lower_limit,upper_limit)+1    
    qi = expm(Q1)[start_index_Q1,obs_index_Q1] # P(X_T=x_T,X_[]<u+1|X_0=x)
    bi = expm(Q2)[start_index_Q2,obs_index_Q2] # P(X_T=x_T|X_0=x)
    di = bi-qi # P(X_T=x_T,X_[]=>u+1|X_0=x)=P(X_T=x_T,X_[]>u|X_0=x)
    return(di)
}

lower = 90
upper = 105
obs = 95
lower_limit = 0
upper_limit = 300
x = c(100)
theta = c(200,200)

df = get_Q(x,theta,obs,lower,upper,lower_limit,upper_limit)
bridge_index = df$bridge_index
target_index = df$target_index
expQ = df$expQ
plot_Q(x,theta,obs,lower,upper,lower_limit,upper_limit)

mu = get_mu("BDI_2",x,theta)
sig = sqrt(get_var("BDI_2",x,theta))

obs = 105
x = 100
get_above_maximum_sigma_wt(upper-x+1,mu,sig,obs-x,1)
get_Q_max(x,theta,obs,lower,upper,lower_limit,upper_limit)

tau = 0.1
data = sim_data("BDI_2",x,theta, tout, lower, upper, tau)$data
data_chain = get_chain(data,tout)
plot_chain(data_chain,c("Population"),c("blue"),ylim=c(50,150))
add_box(lower,upper,tout,tau,"blue")
add_obs(tout,obs,"blue")

Vget_above_maximum_sigma_wt = Vectorize(get_above_maximum_sigma_wt,vectorize.args = 'wt')
exact_brown = function(x,mu,sig,wt,tout){
    integrate(Vget_above_maximum_sigma_wt,wt-1/2,wt+1/2,mu=mu,sig=sig,x=x,tout=tout)$value
}
vexact_brown = Vectorize(exact_brown,vectorize.args = c('x','wt'))

sum(vexact_brown(upper-(lower_limit:upper_limit)+1,mu,sig,obs-(lower_limit:upper_limit),1))

