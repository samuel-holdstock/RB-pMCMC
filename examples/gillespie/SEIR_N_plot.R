N = 1000
S0 = 997
E0 = 0
I0 = 3
R0 = N-S0-I0
beta = 0.4
alpha = 0.25
gamma = 0.04
total_time = 125

data = sim_data(c(S0,E0,I0),c(N,beta,alpha,gamma), SEIR.S, total_time, c(1050,100,20), c(900,0,0), 0.1,"SEIR_N")
data = cbind(data,N - apply(data[,-1],1,sum))
data_chain = get_chain(data,total_time)

plot(data_chain[,c(1,2)],col='blue',type='l')
lines(data_chain[,c(1,3)],col='orange')
lines(data_chain[,c(1,4)],col='darkgreen')
lines(data_chain[,c(1,5)],col='red')


