source("Rlibraries")

sourceCpp("main.cpp")
run_algorithm(c(1000),c(0.1,0.1,50), BDI.S, 1, c(1100), c(900), 0.1,"BDI")
run_algorithm(c(1000,5),c(0.01,0.25), SIR.S, 1, c(1050,20), c(900,0), 0.1,"SIR")
run_algorithm(c(1000,0,5),c(.018,1/12,1/3), SEIR.S, 1, c(1050,20,20), c(900,0,0), 0.1,"SEIR")

plot(sim_data(c(1000),c(0.05,0.1,50), BDI.S, 1, c(1100), c(900), 0.1,"BDI"))

sim_data(c(1000),c(0.05,0.1,50), BDI.S, 1, c(1100), c(900), 0.1,"BDI")
plot(sim_data(c(1000,5),c(0.01,0.25), SIR.S, 1, c(1050,20), c(900,0), 0.1,"SIR")[,-3])
plot(sim_data(c(1000,5),c(0.01,0.25), SIR.S, 1, c(1050,20), c(900,0), 0.1,"SIR")[,-2])
sim_data(c(1000,0,5),c(.018,1/12,1/3), SEIR.S, 1, c(1050,20,20), c(900,0,0), 0.1,"SEIR")

