N = 1000
S0 = 997
I0 = 3
R0 = N-S0-I0
beta = 0.4
gamma = 0.04
tout = 125

x = c(S0,I0)
theta = c(N,beta,gamma)

obs = c(0,13)

lower = c(0,0)
upper = c(2,50)

lower_limit = c(0,0)
upper_limit = c(5,100)
