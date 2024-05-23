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

#endif