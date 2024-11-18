#ifndef SWARM_H
#define SWARM_H

#include <iostream>
#include <random>
#include <memory>
#include <unordered_map>
#include <string>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class Swarm {
  public: 
    Swarm();
    int time;
    std::vector<std::vector<std::pair<std::pair<arma::vec,int>,double>>> swarm;
    std::vector<std::pair<std::pair<arma::vec,int>,double>> pairs;
    std::vector<int> parent_indices;
    // std::vector<std::pair<arma::vec,int>> states;
    // std::vector<double> probabilities;
    // std::vector<arma::vec> path;
    void add_state(const std::pair<arma::vec,int> &state,const double &probability);
    std::pair<arma::vec,int> get_parent(int parent_index);
    void push_state();
    void pop_state();
    arma::mat get_path();
};


#endif
