str = "SIR_N"
N = 100
beta = 0.4
gamma = 0.04
theta = c(N,beta,gamma)
S = SIR.S
tau = 1
lower = c(1,1)
upper = c(20,20)

Q = get_coffin_matrix("SIR_N",lower,upper,theta)
Qs = as(Q,"sparseMatrix")

lower_mat = matrix(c(1,1,
                     1,1,
                     1,1,
                     1,1),ncol=2,byrow=T)
upper_mat = matrix(c(5,5,
                     10,10,
                     20,20,
                     30,30),ncol=2,byrow=T)

obs_index = 4
benchmarkQ = benchmark_exp_time(str,lower_mat,upper_mat,theta,S,obs_index,tau,F)
benchmarkQs = benchmark_exp_time(str,lower_mat,upper_mat,theta,S,obs_index,tau,T)
plot(benchmarkQ$times,main="Time to calculate Exp(Qtau)",xlab='Test case',ylab='Median (micro secs)',type='o',col='blue')
lines(benchmarkQs$times,type='o',col='red')
legend("topleft",legend=c("Dense","Sparse"),col=c("blue","red"),inset=0.05,lty=c(1,1))
