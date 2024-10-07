#ifndef GILLESPIE_H
#define GILLESPIE_H

#include "get_Q.h"
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

Rcpp::List gillespie_alg(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

Rcpp::List gillespie_alg_entire(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

Rcpp::List gillespie_alg_frac(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

Rcpp::List gillespie_alg_entire_frac(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

arma::mat add_row_time(arma::mat data, arma::vec xcurr, double tnext);

#endif