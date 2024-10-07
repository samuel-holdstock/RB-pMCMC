#ifndef MCMC_H
#define MCMC_H

#include "main.h"
#include "rates.h"
#include "brownian.h"

#include <math.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

arma::vec propose_psi(const arma::vec &theta);
double get_likelihood_RB(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &num_particles);
double get_likelihood_frac(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &num_particles);
double frac_list(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &obs_list, int M);
arma::mat mvrnormArma(int n, arma::mat sigma);
double get_prior(const arma::vec log_theta);

#endif