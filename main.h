#ifndef MAIN_H
#define MAIN_H

#include "libraries.h"

#include "utils/rexpQ.h"
#include "get_Q.h"
#include "debug.h"
#include "gillespie.h"
#include "rates.h"
#include "brownian.h"
#include "exact.h"
#include "mcmc.h"
#include "utils/swarm.h"

std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> get_rate_function(const std::string &str);

arma::mat get_S(const std::string &str);

arma::vec get_rate(const std::string &str, const arma::vec &x, const arma::vec &theta);

Rcpp::List run_algorithm(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau);

Rcpp::List run_algorithm_frac(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout);

Rcpp::List sim_data(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau);

Rcpp::List sim_data_frac(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout);

arma::mat get_coffin_matrix(
std::string str,
arma::vec lower,
arma::vec upper,
arma::vec theta);
double get_f(const arma::vec &xt, const arma::vec &yt);

double get_estimate(const Rcpp::List &estimate, const arma::vec &lower, const arma::vec &upper, const arma::mat &P, const arma::vec &obs);

double get_estimate_frac(const Rcpp::List &estimate, const arma::vec &obs);

uint64_t timeSinceEpochMillisec();
uint64_t timeSinceEpochNanosec();

double RB(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau,
                  const arma::vec &obs, const int &M);

double RB_list(const std::string &str, const arma::vec &x0, const arma::vec &theta, const arma::vec &tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, const arma::vec &tau_list,
                  const arma::mat &obs_list, const int &M);
double RB_particle_filter(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M);
Rcpp::List RB_particle_filter_gibbs(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M);

double frac(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &obs, const int &M);
double frac_list(const std::string &str, const arma::vec &x0, const arma::vec &theta, const arma::vec &tout_list,
                  const arma::mat &obs_list, const int &M);
double frac_particle_filter(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &obs_list, int M);
Rcpp::List frac_particle_filter_gibbs(const std::string &str, const arma::vec &x0, const arma::vec &theta, const arma::vec &tout_list,
                  const arma::mat &obs_list, const int &M);

std::unordered_map<std::string,int> state_to_index_map(const arma::mat &i2s_map);
arma::mat index_to_state_map(const arma::vec &lower,const arma::vec &upper);
std::string vec_to_string(const arma::vec &vec);
arma::vec string_to_vec(const std::string &str);
int state_to_index_map_convert(const arma::vec &state, std::unordered_map<std::string,int> &s2i_map);
#endif


