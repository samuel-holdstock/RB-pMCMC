#include <iostream>
#include <random>
#include <memory>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "gillespie.h"
#include "rates.h"
#include "get_Q.h"

std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> initialize_algorithm(std::string str){
  ModelDictionary model_dict;
  auto model = model_dict.get_model(str);
  return [model](const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
    return model->get_rates(x,theta);
  };

//   std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = ->get_rate_function();
//   Rcpp::NumericVector x={10};
//   Rcpp::NumericVector theta={0.1,0.1,0.1};
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
  std::cout<<"A";
  return(gillespie_alg_entire(x0, theta, S, tout, upper, lower, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix get_coffin_matrix(
Rcpp::NumericVector lower,
Rcpp::NumericVector upper,
Rcpp::NumericVector theta,
Rcpp::NumericMatrix S,
std::string str
){
  std::cout<<"A";
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = initialize_algorithm(str);
  std::cout<<"B";
  return(Qmat(lower,upper,theta,S,rates_function));
}
