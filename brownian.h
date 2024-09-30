#ifndef BROWNIAN_H
#define BROWNIAN_H

#include "main.h"
#include "rates.h"

#include <math.h>
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>

double get_above_maximum_wt(double x,double mu, double wt, double tout);
double get_above_maximum_sigma_wt(double x, double mu, double sig, double wt, double tout);
double get_below_minimum_sigma_wt(double x, double mu, double sig, double wt, double tout); 
double get_probability_hit(double start, double target, double tout, double mu, double sig);
Rcpp::NumericVector get_mu(const std::string &model_name, const Rcpp::NumericVector &thetas, const Rcpp::NumericVector &start);
Rcpp::NumericVector get_var(const std::string &model_name, const Rcpp::NumericVector &thetas, const Rcpp::NumericVector &start);
double get_variance_brownian(const std::string &model_name, const Rcpp::NumericVector &thetas, const Rcpp::NumericVector &start);
double get_var_big_box_tau(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &target, const double &tout, const double tau);
double get_PVR_big_box_tau(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &target, const double &tout, const double tau);
Rcpp::List get_box_brownian_fast_list(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &tout_list, const Rcpp::NumericVector &tau_list, const Rcpp::NumericVector &start, const Rcpp::NumericMatrix &target_list,
const double &pvr_goal);
double get_tau(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &target, const double &tout, const double &goal);
Rcpp::NumericVector get_tau_list(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const double &goal);
Rcpp::List get_box_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target,
const double &pvr_goal);

#endif