#ifndef GET_Q_H
#define GET_Q_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

void calculate_Q(Rcpp::NumericVector upper, Rcpp::NumericVector lower);

//[[Rcpp::export]]
int state_to_index(
Rcpp::NumericVector state,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper);

//[[Rcpp::export]]
Rcpp::NumericVector index_to_state(
int index,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper);

Rcpp::NumericMatrix Qmat(
Rcpp::NumericVector lower,
Rcpp::NumericVector upper,
Rcpp::NumericVector theta,
Rcpp::NumericMatrix S,
std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function
);

bool inBox(
Rcpp::NumericVector state,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper
);

#endif