#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <random>
#include <memory>

#include "utils/rexpQ.h"
#include "get_Q.h"
#include "gillespie.h"
#include "rates.h"
#include "brownian.h"
#include "exact.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> get_rate_function(std::string str);

Rcpp::NumericMatrix get_S(std::string str);

Rcpp::NumericVector get_rate(std::string str, const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);

Rcpp::List run_algorithm(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau);

Rcpp::List run_algorithm_frac(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout);

Rcpp::List sim_data(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau);

Rcpp::List sim_data_frac(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout);

Rcpp::NumericMatrix get_coffin_matrix(
std::string str,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper,
Rcpp::NumericVector theta);

double get_estimate(const Rcpp::List &estimate, const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, const arma::mat &P, const Rcpp::NumericVector &obs);

double get_estimate_frac(const Rcpp::List &estimate, const Rcpp::NumericVector &obs);

double RB(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau,
                  const Rcpp::NumericVector &obs, int M);

double frac(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &obs, int M);

#endif


