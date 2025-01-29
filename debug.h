#ifndef DEBUG_H
#define DEBUG_H

#include "main.h"


class DebugRB {
  public: 
  arma::mat x;
  double log_prob;
  arma::vec inS;
  arma::vec hitS;
  arma::mat xt_data;
  arma::mat xttau_data;
  arma::mat presample;
  arma::vec weights;
  arma::vec particle_weights;
  uint64_t exp_time;
  uint64_t entire_time;
  uint64_t gillespie_time;
  uint64_t weight_time;
  uint64_t resample_time;
  uint64_t enter_box_check_time = 0.0;
  uint64_t reaction_time = 0.0;
  uint64_t left_box_check_time = 0.0;
  uint64_t no_reaction_check_time = 0.0;
  uint64_t init_time=0.0;
  Rcpp::List get_tuple();
};
class DebugFrac {
  public: 
  arma::mat x;
  double log_prob;
  arma::mat xt_data;
  arma::mat presample;
  arma::vec weights;
  uint64_t entire_time;
  uint64_t gillespie_time;
  uint64_t weight_time;
  uint64_t resample_time;
  Rcpp::List get_tuple();
};

#endif