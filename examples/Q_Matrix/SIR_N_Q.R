str = "SIR_N"
N = 100
beta = 0.4
gamma = 0.04
theta = c(N,beta,gamma)
S = SIR.S
tau = 1
lower = c(1,1)
upper = c(20,20)

Q = get_coffin_matrix(lower,upper,theta,S,"SIR_N")
Qs = as(Q,"sparseMatrix")

lower_mat = matrix(c(1,1,
                     1,1,
                     1,1,
                     1,1),ncol=2,byrow=T)
upper_mat = matrix(c(5,5,
                     10,10,
                     20,20,
                     30,30),ncol=2,byrow=T)

benchmark = benchmark_exp_time(lower_mat,upper_mat,theta,S,"SIR_N",4,tau)
plot(benchmark$times,xlab='Upper',ylab='Median time (micro)',type='o')
