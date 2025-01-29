#ifndef MCMC_H
#define MCMC_H

#include "main.h"

arma::vec propose_psi(const arma::vec &theta);
double get_likelihood_RB(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::mat &lower_list, const arma::mat &upper_list, const arma::vec &tout_list, const int &num_particles);
double get_likelihood_frac(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &num_particles);
Rcpp::List get_likelihood_RB_gibbs(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::mat &lower_list, const arma::mat &upper_list, const arma::vec &tout_list, const int &num_particles);
Rcpp::List get_likelihood_frac_gibbs(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &num_particles);
arma::mat mvrnormArma(int n, arma::mat sigma);
double get_prior(const arma::vec &transformed_theta);
arma::vec transform(const arma::vec &x);
arma::vec inv_transform(const arma::vec &x);
#endif