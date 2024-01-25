#ifndef GILLESPIE_H
#define GILLESPIE_H

#include "get_Q.h"

Rcpp::List gillespie_alg(const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, const Rcpp::NumericMatrix &S, double tout,
                            const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau,
                            std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function);

Rcpp::List gillespie_alg_entire(const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, const Rcpp::NumericMatrix &S, double tout,
                            const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau,
                            std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function);

Rcpp::NumericMatrix add_row_time(Rcpp::NumericMatrix data, Rcpp::NumericVector xcurr, double tnext);

#endif