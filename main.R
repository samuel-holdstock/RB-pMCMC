source("utils/Rlibraries.R")

sourceCpp("main.cpp")

source("examples/gillespie/SIR_N_plot.R")
source("examples/gillespie/SEIR_N_plot.R")

lower = c(10)
upper = c(20)

theta = c(0.1,0.1,0.1)
S = BDI.S

test(upper,lower)
state_to_index(c(1,0,2),c(1,0,2),c(3,5,3))
index_to_state(34,c(1,0,2),c(3,5,3))

m = test(lower,upper,theta,S,"BDI")
test(lower,upper,theta,S,"BDI")
m
expm(m)
SS_exp_Q(m,prec=1e-15)
