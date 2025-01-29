#ifndef BROWNIAN_H
#define BROWNIAN_H

#include "main.h"

double get_above_maximum_wt(double x,double mu, double wt, double tout);
double get_above_maximum_sigma_wt(double x, double mu, double sig, double wt, double tout);
double get_below_minimum_sigma_wt(double x, double mu, double sig, double wt, double tout); 
double get_probability_hit(double start, double target, double tout, double mu, double sig);
arma::vec get_mu(const std::string &model_name, const arma::vec &thetas, const arma::vec &start);
arma::vec get_var(const std::string &model_name, const arma::vec &thetas, const arma::vec &start);
arma::mat get_covar(const std::string &model_name, const arma::vec &thetas, const arma::vec &start);
double get_variance_brownian(const std::string &model_name, const arma::vec &thetas, const arma::vec &start);
double get_var_big_box_tau(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::vec &target, const double &tout, const double tau);
double get_PVR_big_box_tau(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::vec &target, const double &tout, const double tau);
double get_VRF_big_box_tau(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::vec &target, const double &tout, const double tau);
Rcpp::List get_box_brownian_fast_list(const std::string &model_name, const arma::vec &thetas, 
const arma::vec &tout_list, const arma::vec &tau_list, const arma::vec &start, const arma::mat &target_list,
const double &beta);
double get_tau(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::vec &target, const double &tout, const double &alpha);
arma::vec get_tau_list(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::mat &obs_list, const arma::vec &tout_list, const double &alpha);
Rcpp::List get_box_brownian_fast(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const double &tau, const arma::vec &start, const arma::vec &target,
const double &beta);

#endif