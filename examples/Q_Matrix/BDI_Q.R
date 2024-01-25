str = "BDI"
lambda = 0.1
mu = 0.1
gamma = 0.1
theta = c(lambda,mu,gamma)
S = BDI.S
tau = 1
lower = c(1)
upper = c(100)

Q = get_coffin_matrix("BDI",lower,upper,theta)
Qs = as(Q,"sparseMatrix")

#index_to_state(c(1),lower,upper)
#state_to_index(x,lower,upper)

#get_rate("BDI",100,theta)
#get_S("BDI")

lower_mat = matrix(rep(1,20),ncol=1,byrow=T)
upper_mat = matrix(seq(from=10,by=10,length.out=20),ncol=1,byrow=T)

obs_index = 4
benchmarkQ = benchmark_exp_time(str,lower_mat,upper_mat,theta,S,obs_index,tau,F)
benchmarkQs = benchmark_exp_time(str,lower_mat,upper_mat,theta,S,obs_index,tau,T)
plot(benchmarkQ$times,main="Time to calculate Exp(Qtau)",xlab='Test case',ylab='Median (micro secs)',type='o',col='blue')
lines(benchmarkQs$times,type='o',col='red')
legend("topleft",legend=c("Dense","Sparse"),col=c("blue","red"),inset=0.05,lty=c(1,1))

