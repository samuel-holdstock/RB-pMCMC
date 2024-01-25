#include <iostream>
#include <random>
#include <memory>

#include "utils/rexpQ.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "get_Q.h"
#include "gillespie.h"
#include "rates.h"

std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> get_rate_function(std::string str){
  ModelDictionary model_dict;
  auto model = model_dict.get_model(str);
  return [model](const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
    return model->get_rates(x,theta);
  };
}
//[[Rcpp::export]]
Rcpp::NumericMatrix get_S(std::string str){
  ModelDictionary model_dict;
  auto model = model_dict.get_model(str);
  return(model->S);
}
//[[Rcpp::export]]
Rcpp::NumericVector get_rate(std::string str, const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  return(get_rate_function(str)(x,theta));
}


// PREDEFINED MODELS
// [[Rcpp::export]]
Rcpp::List run_algorithm(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  return(gillespie_alg(x0, theta, S, tout, lower, upper, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::List sim_data(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  return(gillespie_alg_entire(x0, theta, S, tout, lower, upper, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::NumericMatrix get_coffin_matrix(
std::string str,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper,
Rcpp::NumericVector theta){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  return(Qmat(lower,upper,theta,S,rates_function));
}

double get_estimate(const Rcpp::List &estimate, const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, const arma::mat &P, const Rcpp::NumericVector &obs){  
  int hitObs = 1;
  Rcpp::NumericVector xt_data = estimate["xt_data"];
  Rcpp::NumericVector xttau_data = estimate["xttau_data"];
  bool inS = estimate["inS"]; 
  for(int i=0;i<xt_data.length();++i){
    hitObs *= (xt_data[i]==obs[i]);
  }
  return(hitObs*inS + P[state_to_index(xttau_data,lower,upper)]);
}

//[[Rcpp::export]]
double RB(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau,
                  const Rcpp::NumericVector &obs, int M){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  Rcpp::NumericMatrix Q = Rcpp::transpose(get_coffin_matrix(str,lower,upper,theta));
  arma::mat v(Q.nrow(),1);    
  int obs_index = state_to_index(obs,lower,upper);
  v[obs_index] = 1;
  arma::mat P = vT_exp_Q(v,Q*tau,1e-15,false,true,false);
  Rcpp::List estimate;
  double estimator = 0;
  for(int i=0;i<M;++i){
    estimate = run_algorithm(str,x0,theta, tout, lower, upper, tau);
    estimator += get_estimate(estimate,lower,upper,P,obs);
  }
  return(estimator/M);
}



