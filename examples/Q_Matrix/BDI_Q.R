str = "BDI"
lambda = 0.1
mu = 0.1
gamma = 0.1
theta = c(lambda,mu,gamma)
S = BDI.S
tau = 1
lower = c(1)
upper = c(100)

Q = get_coffin_matrix(lower,upper,theta,S,"BDI")
Qs = as(Q,"sparseMatrix")

x = index_to_state(c(1),lower,upper)
BDI_get_rates(x,theta)

lower_mat = matrix(rep(1,10),ncol=1,byrow=T)
upper_mat = matrix(seq(from=10,by=10,length.out=10),ncol=1,byrow=T)

obs_index = 4
benchmark = benchmark_exp_time(lower_mat,upper_mat,theta,S,str,obs_index,tau)
plot(benchmark$times,xlab='Upper',ylab='Median time (micro)',type='o')




