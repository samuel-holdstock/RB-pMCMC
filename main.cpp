#include "main.h"

std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> get_rate_function(std::string str){
  auto model = model_dict.get_model(str);
  return [model](const arma::vec &x, const arma::vec &theta){
    return model->get_rates(x,theta);
  };
}
//[[Rcpp::export]]
arma::mat get_S(std::string str){
  auto model = model_dict.get_model(str);
  return(model->S);
}
//[[Rcpp::export]]
arma::vec get_rate(std::string str, const arma::vec &x, const arma::vec &theta){
  return(get_rate_function(str)(x,theta));
}


// PREDEFINED MODELS
// [[Rcpp::export]]
Rcpp::List run_algorithm(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg(x0, theta, S, tout, lower, upper, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::List run_algorithm_frac(std::string str, const arma::vec &x0, const arma::vec &theta, double tout){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_frac(x0, theta, S, tout, rates_function));
}

// [[Rcpp::export]]
Rcpp::List sim_data(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_entire(x0, theta, S, tout, lower, upper, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::List sim_data_frac(std::string str, const arma::vec &x0, const arma::vec &theta, double tout){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_entire_frac(x0, theta, S, tout, rates_function));
}

// [[Rcpp::export]]
arma::mat get_coffin_matrix(
std::string str,
arma::vec lower,
arma::vec upper,
arma::vec theta){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(Qmat(lower,upper,theta,S,rates_function));
}

//[[Rcpp::export]]
double get_estimate(const Rcpp::List &estimate, const arma::vec &lower, const arma::vec &upper, const arma::mat &P, const arma::vec &obs){  
  int hitObs = 1;
  arma::vec xt_data = estimate["xt_data"];
  arma::vec xttau_data = estimate["xttau_data"];
  bool inS = estimate["inS"]; 
  for(int i=0;i<xt_data.n_elem;++i){
    hitObs *= (xt_data[i]==obs[i]);
  }
  //std::cout<<"Hit xt:"<<hitObs<<", Not inS:"<<(1-inS)<<", Q:"<<P[state_to_index(xttau_data,lower,upper)]<<std::endl;
  double est = hitObs*(1-inS) + P[state_to_index(xttau_data,lower,upper)];
  return(est);
}

//[[Rcpp::export]]
double get_estimate_frac(const Rcpp::List &estimate, const arma::vec &obs){  
  int hitObs = 1;
  arma::vec xt_data = estimate["xt_data"];
  for(int i=0;i<xt_data.n_elem;++i){
    hitObs *= (xt_data[i]==obs[i]);
  }
  double est = hitObs;
  //std::cout<<est<<std::endl;
  return(est);
}

//[[Rcpp::export]]
double RB(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau,
                  const arma::vec &obs, int M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  arma::mat Q = get_coffin_matrix(str,lower,upper,theta).t();
  arma::mat v(Q.n_rows,1);    
  int obs_index = state_to_index(obs,lower,upper);
  v[obs_index] = 1;
  arma::mat P = vT_exp_Q(v,Rcpp::wrap(Q*tau),1e-20,false,true,false);
  Rcpp::List estimate;
  double estimator = 0;
  for(int i=0;i<M;++i){
    estimate = run_algorithm(str,x0,theta, tout, lower, upper, tau);
    estimator += get_estimate(estimate,lower,upper,P,obs);
  }
  return(estimator/M);
}

//[[Rcpp::export]]
double RB_list(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec lower(num_species);
  arma::vec upper(num_species);
  arma::vec obs(num_species);
  double tout;
  double tau;
  arma::vec x = x0; // copy
  double prev_time = 0;
  double log_prob = 0;
  for(int i=0; i<num_obs; ++i){
    lower = lower_list.row(i);
    upper = upper_list.row(i);
    obs = obs_list.row(i);
    tout = tout_list(i)-prev_time;
    tau = tau_list[i];
    log_prob += log(RB(str,x,theta,tout,lower,upper,tau,obs,M));
    if(isinf(log_prob)){
      return(log_prob);
    }
    x = obs; // copy
    prev_time = tout_list(i);
  }
  return(log_prob);
}

//[[Rcpp::export]]
double frac(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &obs, int M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  Rcpp::List estimate;
  double estimator = 0;
  for(int i=0;i<M;++i){
    estimate = run_algorithm_frac(str,x0,theta, tout);
    estimator += get_estimate_frac(estimate,obs);
  }
  return(estimator/M);
}

//[[Rcpp::export]]
double frac_list(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &obs_list, int M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec obs(num_species);
  double tout;
  arma::vec x = x0; // copy
  double prev_time = 0;
  double log_prob = 0;
  for(int i=0; i<num_obs; ++i){
    obs = obs_list.row(i);
    tout = tout_list(i)-prev_time;
    log_prob += log(frac(str,x,theta,tout,obs,M));
    if(isinf(log_prob)){
      return(log_prob);
    }
    x = obs; // copy
    prev_time = tout_list(i);
  }
  return(log_prob);
}


