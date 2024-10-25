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

// [[Rcpp::export]]
double get_f(const arma::vec &xt, const arma::vec &yt){
  // Inexact Uniform:
  // int box_halfwidth = 5;
  // if(arma::sum(arma::abs(xt-yt))<=box_halfwidth){
  //   return(1/pow(2*box_halfwidth+1,xt.size()));
  // }

  // Inexact Missing Prey:  
  // if(fabs(xt(1)-yt(0))<0.5){
  //   return(1);
  // }
  // return(0);

  // Inexact Gaussian Noise
  arma::vec zt = xt-yt;
  int d = zt.size();
  double f = 1;
  double noise_sd = 2;
  for(int i=0;i<d;++i){
    f = f*R::dnorm(zt(i),0,noise_sd,0);
  }
  return(f);
}

// [[Rcpp::export]]
arma::vec get_ft(const arma::vec &lower, const arma::vec &upper, const arma::vec &yt, const arma::mat &P){
  int num_param = lower.n_elem;
  int total_points = P.n_rows;
  arma::vec f(total_points);
  for(int obs_index=0;obs_index<total_points-1;++obs_index){
    f(obs_index) = get_f(index_to_state(obs_index,lower,upper),yt);
  }
  f(total_points-1)=0;
  return(f);
}


// [[Rcpp::export]]
arma::vec get_fP(const arma::vec &lower, const arma::vec &upper, const arma::vec &yt, const arma::mat &P){
  arma::vec f = get_ft(lower,upper,yt,P);
  return((f.t()*P).t());
}

bool obs_in_box(const arma::vec &lower, const arma::vec &upper, const arma::vec &obs, const arma::mat &P){
  int index = state_to_index(obs,lower,upper);
  if(index==P.n_rows-1){
    std::cout<<"WARNING: OBS NOT IN BOX"<<std::endl;
    return(false);
  }
  return(true);
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
  // std::cout<<"Hit xt:"<<hitObs<<", Not inS:"<<(1-inS)<<", Q:"<<P[state_to_index(xttau_data,lower,upper)]<<std::endl;
  // std::cout<<"XT: "<<xt_data<<std::endl;
  // std::cout<<"X_T-TAU: "<<xttau_data<<std::endl;
  // std::cout<<"P: "<<P<<std::endl;
  double est = hitObs*(1-inS) + P[state_to_index(xttau_data,lower,upper)];
  return(est);
}

//[[Rcpp::export]]
double get_estimate_noise_P(const Rcpp::List &estimate, const arma::vec &lower, const arma::vec &upper, const arma::mat &P, const arma::vec &obs, const arma::vec &fP){  
  arma::vec xt_data = estimate["xt_data"];
  arma::vec xttau_data = estimate["xttau_data"];
  bool inS = estimate["inS"]; 
  int tau_index = state_to_index(xttau_data,lower,upper);
  // std::cout<<"X_{T-tau}: "<<xttau_data;
  // std::cout<<"X_{T}: "<<xt_data;
  // std::cout<<"inS: "<<inS<<std::endl;
  // std::cout<<"Index: "<<tau_index<<std::endl;
  // std::cout<<"fP: "<<fP;
  // std::cout<<"f(xt): "<<get_f(xt_data,obs)<<std::endl;
  // std::cout<<"Hit xt:"<<hitObs<<", Not inS:"<<(1-inS)<<", Q:"<<P[state_to_index(xttau_data,lower,upper)]<<std::endl;
  double est = get_f(xt_data,obs)*(1-inS) + fP(tau_index);
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
double get_estimate_frac_noise(const Rcpp::List &estimate, const arma::vec &obs){  
  arma::vec xt_data = estimate["xt_data"];
  return(get_f(xt_data,obs));
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
  if(!obs_in_box(lower,upper,obs,Q)){
    return(0);
  }
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
std::string vec_to_string(arma::vec vec){
  std::ostringstream oss;
  oss.precision(5);
  for(int i=0; i<vec.n_elem;++i){
    oss<<vec[i]<<",";
  }
  return(oss.str());
}

void add_to_map(std::string key, double weight, std::unordered_map<std::string,double> &map){
  if(map.find(key)==map.end()){
    map[key] = 0;
  }
  map[key] += weight;
}

//[[Rcpp::export]]
arma::vec string_to_vec(const std::string &str){
  std::vector<double> values;
  std::istringstream iss(str);
  std::string token;
  while(std::getline(iss,token,',')){
    if(!token.empty()){
      values.push_back(std::stod(token));
    }
  }
  return(arma::vec(values));
}

arma::mat sample_x(const arma::vec &resample_vec_xttau, std::unordered_map<std::string,double> &resample_vec_xt, const arma::vec &lower, const arma::vec &upper, const int &M){
  int xttau_size = resample_vec_xttau.size(); 
  int xt_size = resample_vec_xt.size(); 
  double normalise_constant = 0;
  int num_species = lower.size();

  for(int i=0; i<resample_vec_xttau.size()-1;++i){
    add_to_map(vec_to_string(index_to_state(i,lower,upper)),resample_vec_xttau(i),resample_vec_xt);
  }
  for(const auto &entry : resample_vec_xt){
    normalise_constant += entry.second;
    // std::cout<<"Normalise "<<entry.first<<": "<<normalise_constant<<std::endl;
  }
  std::vector<double> probabilities;
  std::vector<string> keys;
  
  for(auto &entry : resample_vec_xt){
    entry.second /= normalise_constant;
    if(probabilities.empty()){
      probabilities.push_back(entry.second);
    }
    else{
      probabilities.push_back(entry.second+probabilities.back());
    }
    keys.push_back(entry.first);
    // std::cout<<"Key: "<<entry.first<<" Value: "<<entry.second<<std::endl;
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0,1.0);
  // arma::mat samples(K,num_species);
  std::vector<string> samples;
  samples.reserve(M);
  double rand_num;
  for(int i=0;i<M;++i){
    rand_num = dis(gen);
    auto it = std::lower_bound(probabilities.begin(),probabilities.end(),rand_num);
    int index = std::distance(probabilities.begin(),it); 
    samples.push_back(keys[index]);
    // std::cout<<"Sample "<<i<<": "<<keys[index]<<std::endl;
  }  
  arma::mat next_x(M,num_species);
  for(int i=0;i<M;++i){
    next_x.row(i) = string_to_vec(samples[i]).t();
  }
  return(next_x);
}

std::pair<arma::mat,double> RB_noise(std::string str, const arma::mat &x0, const arma::vec &theta, double tout,
                  const arma::vec &lower, const arma::vec &upper, double tau,
                  const arma::vec &obs, int M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  arma::mat Q = get_coffin_matrix(str,lower,upper,theta).t();
  arma::mat P = arma::expmat(Q*tau);
  arma::vec f = get_ft(lower,upper,obs,P);
  int tau_index;
  Rcpp::List estimate;
  std::unordered_map<std::string,double> resample_vec_xt;
  arma::vec resample_vec_xttau(f.size());
  double estimator = 0;
  bool inS;
  arma::vec x;
  for(int i=0;i<M;++i){
    x = x0.row(i).t();
    estimate = run_algorithm(str,x,theta, tout, lower, upper, tau);
    estimator += get_estimate_noise_P(estimate,lower,upper,P,obs,(f.t()*P).t());
    inS = estimate["inS"];
    add_to_map(vec_to_string(estimate["xt_data"]),get_f(estimate["xt_data"],obs)*(1-inS),resample_vec_xt);
    resample_vec_xttau = resample_vec_xttau + f%P.col(state_to_index(estimate["xttau_data"],lower,upper));
  }
  arma::mat next_x = sample_x(resample_vec_xttau,resample_vec_xt,lower,upper,M);
  std::pair<arma::mat,double> result = std::make_pair(next_x, estimator/M);
  return(result);
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
double RB_particle_filter(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec lower(num_species);
  arma::vec upper(num_species);
  arma::vec obs(num_species);
  double tout;
  double tau;
  arma::mat x = arma::repmat(x0.t(),M,1);
  double prev_time = 0;
  double log_prob = 0;
  std::pair<arma::mat,double> result;
  for(int i=0; i<num_obs; ++i){
    lower = lower_list.row(i).t();
    upper = upper_list.row(i).t();
    obs = obs_list.row(i).t();
    tout = tout_list(i)-prev_time;
    tau = tau_list[i];
    result = RB_noise(str,x,theta,tout,lower,upper,tau,obs,M);
    log_prob += log(result.second);
    if(isinf(log_prob)){
      return(log_prob);
    }
    x = result.first; // copy
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
double frac_noise(std::string str, const arma::vec &x0, const arma::vec &theta, double tout,
                  const arma::vec &obs, int M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  Rcpp::List estimate;
  double estimator = 0;
  for(int i=0;i<M;++i){
    estimate = run_algorithm_frac(str,x0,theta, tout);
    estimator += get_estimate_frac_noise(estimate,obs);
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

//[[Rcpp::export]]
double frac_particle_filter(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
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
    log_prob += log(frac_noise(str,x,theta,tout,obs,M));
    if(isinf(log_prob)){
      return(log_prob);
    }
    x = obs; // copy
    prev_time = tout_list(i);
  }
  return(log_prob);
}
