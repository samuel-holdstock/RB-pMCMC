#ifndef MCMC_H
#define MCMC_H

#include "main.h"
#include "rates.h"
#include "brownian.h"

#include <math.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

Rcpp::NumericVector propose_psi(const Rcpp::NumericVector &theta);
double get_likelihood_RB(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &theta, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const int &num_particles);
double get_likelihood_frac(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &theta, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const int &num_particles);
double frac_list(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, Rcpp::NumericVector tout_list,
                  const Rcpp::NumericMatrix &obs_list, int M);
arma::mat mvrnormArma(int n, arma::mat sigma);
double get_prior(const Rcpp::NumericVector log_theta);

#endif