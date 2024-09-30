library("Rcpp")
library("RcppArmadillo")
library("expQ")
library("expm")
library("microbenchmark")
library("scales") 

source("simulation-study/simulation-utils.R")
source("utils/Rfunctions.R")
source("examples/Q_Matrix/Q_Matrix_utils.R")
source("get_Q.r")

plot(z1$log_posterior[,3],type='l')
plot(z1$log_posterior[,3],type='l')

effectiveSize(z1$log_posterior)
