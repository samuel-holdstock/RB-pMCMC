#include <iostream>
#include <random>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "gillespie.h"
#include "rates.h"

std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> initialize_algorithm(std::string str){
    RateDictionary rate_dictionary;
    //rate_dictionary.get_rate_function(str);
    initialize_rates(rate_dictionary);
    return(rate_dictionary.get_rate_function(str));
}

// [[Rcpp::export]]
Rcpp::NumericVector run_algorithm(const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, const Rcpp::NumericMatrix &S, double tout,
                  const Rcpp::NumericVector &upper, const Rcpp::NumericVector &lower, double tau, std::string str){
    std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = initialize_algorithm(str);
    return(gillespie_alg(x0, theta, S, tout, upper, lower, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix sim_data(const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, const Rcpp::NumericMatrix &S, double tout,
                  const Rcpp::NumericVector &upper, const Rcpp::NumericVector &lower, double tau, std::string str){
    std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = initialize_algorithm(str);
    return(gillespie_alg_entire(x0, theta, S, tout, upper, lower, tau, rates_function));
}

