#ifndef EXACT_H
#define EXACT_H

#include "main.h"
#include "get_Q.h"

//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

double get_variance_exact(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper,
const Rcpp::NumericVector &lower_limit, const Rcpp::NumericVector &upper_limit);

Rcpp::NumericMatrix get_contour_exact(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const Rcpp::NumericVector &taus, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list, 
const Rcpp::NumericVector &lower_limit, const Rcpp::NumericVector &upper_limit);

#endif