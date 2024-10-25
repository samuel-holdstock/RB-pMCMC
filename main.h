#ifndef MAIN_H
#define MAIN_H

#include <iostream>
#include <random>
#include <memory>
#include <unordered_map>
#include <string>

#include "utils/rexpQ.h"
#include "get_Q.h"
#include "gillespie.h"
#include "rates.h"
#include "brownian.h"
#include "exact.h"
#include "mcmc.h"

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> get_rate_function(std::string str);

arma::mat get_S(std::string str);

arma::vec get_rate(std::string str, const arma::vec &x, const arma::vec &theta);

Rcpp::List run_algorithm(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau);

Rcpp::List run_algorithm_frac(std::string str, const arma::vec &x0, const arma::vec &theta, double tout);

Rcpp::List sim_data(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau);

Rcpp::List sim_data_frac(std::string str, const arma::vec &x0, const arma::vec &theta, double tout);

arma::mat get_coffin_matrix(
std::string str,
arma::vec lower,
arma::vec upper,
arma::vec theta);

double get_estimate(const Rcpp::List &estimate, const arma::vec &lower, const arma::vec &upper, const arma::mat &P, const arma::vec &obs);

double get_estimate_frac(const Rcpp::List &estimate, const arma::vec &obs);

double RB(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau,
                  const arma::vec &obs, int M);

double RB_list(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M);

double frac(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &obs, int M);

#endif


