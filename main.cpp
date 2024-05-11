#include "main.h"

std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> get_rate_function(std::string str){
  auto model = model_dict.get_model(str);
  return [model](const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
    return model->get_rates(x,theta);
  };
}
//[[Rcpp::export]]
Rcpp::NumericMatrix get_S(std::string str){
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
Rcpp::List run_algorithm_frac(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  return(gillespie_alg_frac(x0, theta, S, tout, rates_function));
}

// [[Rcpp::export]]
Rcpp::List sim_data(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  return(gillespie_alg_entire(x0, theta, S, tout, lower, upper, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::List sim_data_frac(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  return(gillespie_alg_entire_frac(x0, theta, S, tout, rates_function));
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

//[[Rcpp::export]]
double get_estimate(const Rcpp::List &estimate, const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, const arma::mat &P, const Rcpp::NumericVector &obs){  
  int hitObs = 1;
  Rcpp::NumericVector xt_data = estimate["xt_data"];
  Rcpp::NumericVector xttau_data = estimate["xttau_data"];
  bool inS = estimate["inS"]; 
  for(int i=0;i<xt_data.length();++i){
    hitObs *= (xt_data[i]==obs[i]);
  }
  //std::cout<<"Hit xt:"<<hitObs<<", Not inS:"<<(1-inS)<<", Q:"<<P[state_to_index(xttau_data,lower,upper)]<<std::endl;
  double est = hitObs*(1-inS) + P[state_to_index(xttau_data,lower,upper)];
  return(est);
}

//[[Rcpp::export]]
double get_estimate_frac(const Rcpp::List &estimate, const Rcpp::NumericVector &obs){  
  int hitObs = 1;
  Rcpp::NumericVector xt_data = estimate["xt_data"];
  for(int i=0;i<xt_data.length();++i){
    hitObs *= (xt_data[i]==obs[i]);
  }
  double est = hitObs;
  //std::cout<<est<<std::endl;
  return(est);
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
  arma::mat P = vT_exp_Q(v,Q*tau,1e-20,false,true,false);
  Rcpp::List estimate;
  double estimator = 0;
  for(int i=0;i<M;++i){
    estimate = run_algorithm(str,x0,theta, tout, lower, upper, tau);
    estimator += get_estimate(estimate,lower,upper,P,obs);
  }
  return(estimator/M);
}

//[[Rcpp::export]]
double RB_list(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, Rcpp::NumericVector tout_list,
                  const Rcpp::NumericMatrix &lower_list, const Rcpp::NumericMatrix &upper_list, double tau,
                  const Rcpp::NumericMatrix &obs_list, int M){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  int num_obs = tout_list.size();
  Rcpp::NumericVector lower;
  Rcpp::NumericVector upper;
  Rcpp::NumericVector obs;
  double tout;
  Rcpp::NumericMatrix Q;
  arma::mat P;
  Rcpp::List estimate;
  Rcpp::NumericVector estimator(num_obs);
  Rcpp::NumericVector x = clone(x0);
  double prev_time = 0;
  for(int i=0; i<num_obs; ++i){
    lower = lower_list(i,Rcpp::_);
    upper = upper_list(i,Rcpp::_);
    obs = obs_list(i,Rcpp::_);
    tout = tout_list(i)-prev_time;
    Q = Rcpp::transpose(get_coffin_matrix(str,lower,upper,theta));
    arma::mat v(Q.nrow(),1);    
    int obs_index = state_to_index(obs,lower,upper);
    v[obs_index] = 1;
    P = vT_exp_Q(v,Q*tau,1e-20,false,true,false);
    for(int j=0;j<M;++j){
      estimate = run_algorithm(str,x,theta, tout, lower, upper, tau);
      estimator[i] += get_estimate(estimate,lower,upper,P,obs);
    }
    x = clone(obs);
    prev_time = tout;
  }
  double prob = 1;
  for(int i=0; i<num_obs; ++i){
    prob *= estimator[i]/M;
  }
  return(prob);
}

//[[Rcpp::export]]
double frac(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, double tout,
                  const Rcpp::NumericVector &obs, int M){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  Rcpp::List estimate;
  double estimator = 0;
  for(int i=0;i<M;++i){
    estimate = run_algorithm_frac(str,x0,theta, tout);
    estimator += get_estimate_frac(estimate,obs);
  }
  return(estimator/M);
}

//[[Rcpp::export]]
double frac_list(std::string str, const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, Rcpp::NumericVector tout_list,
                  const Rcpp::NumericMatrix &obs_list, int M){
  std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function = get_rate_function(str);
  Rcpp::NumericMatrix S = get_S(str);
  Rcpp::List estimate;
  int num_obs = tout_list.size();
  Rcpp::NumericVector estimator(num_obs);
  Rcpp::NumericVector x = clone(x0);
  Rcpp::NumericVector obs;
  double prev_time = 0;
  for(int i=0;i<num_obs;++i){  
    obs = obs_list(i,Rcpp::_);
    for(int j=0;j<M;++j){
      estimate = run_algorithm_frac(str,x,theta, tout_list(i)-prev_time);
      estimator[i] += get_estimate_frac(estimate,obs);
    }
    x = clone(obs);
    prev_time = tout_list(i);
  }
  double prob = 1;
  for(int i=0;i<num_obs;++i){
    prob *= estimator[i]/M;
  }
  return(prob);
}

