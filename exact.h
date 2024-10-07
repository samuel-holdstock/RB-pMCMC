#ifndef EXACT_H
#define EXACT_H

#include "main.h"
#include "get_Q.h"

//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

double get_variance_exact(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const double &tau, const arma::vec &start, const arma::vec &target, 
const arma::vec &lower, const arma::vec &upper,
const arma::vec &lower_limit, const arma::vec &upper_limit);

arma::mat get_contour_exact(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const arma::vec &taus, const arma::vec &start, const arma::vec &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list, 
const arma::vec &lower_limit, const arma::vec &upper_limit);

#endif