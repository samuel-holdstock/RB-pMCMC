#ifndef GET_Q_H
#define GET_Q_H
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

void calculate_Q(arma::vec upper, arma::vec lower);

//[[Rcpp::export]]
int state_to_index(
arma::vec state,
arma::vec lower,
arma::vec upper);

//[[Rcpp::export]]
arma::vec index_to_state(
int index,
arma::vec lower,
arma::vec upper);

arma::mat Qmat(
arma::vec lower,
arma::vec upper,
arma::vec theta,
arma::mat S,
std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function
);

bool inBox(
arma::vec state,
arma::vec lower,
arma::vec upper
);

#endif