pred = 200
prey = 200

alpha = 0.5
beta = 0.0025
gamma = 0.3
tout = 20

#alpha/beta
#gamma/beta

x = c(pred,prey)
theta = c(alpha,beta,gamma)

obs = c(292,112)
lower = c(250,70)
upper = c(330,150)

lower_limit = c(0,0)
upper_limit = c(30,30)

tout = 1
obs = c(240,190)
lower = c(obs[1]-30,obs[2]-30)
upper = c(obs[1]+30,obs[2]+30)

