N = 1000
S0 = 997
E0 = 0
I0 = 3
R0 = N-S0-I0
beta = 0.4
alpha = 0.25
gamma = 0.04
tout = 125

x = c(S0,E0,I0)
theta = c(N,beta,alpha,gamma)

obs = c(0,0,30)
