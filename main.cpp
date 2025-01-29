#include "main.h"
// #include "libraries.h"
#include <RcppArmadilloExtensions/sample.h>

std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> get_rate_function(const std::string &str){
  auto model = model_dict.get_model(str);
  return [model](const arma::vec &x, const arma::vec &theta){
    return model->get_rates(x,theta);
  };
}
//[[Rcpp::export]]
arma::mat get_S(const std::string &str){
  auto model = model_dict.get_model(str);
  return(model->S);
}
//[[Rcpp::export]]
arma::vec get_rate(const std::string &str, const arma::vec &x, const arma::vec &theta){
  return(get_rate_function(str)(x,theta));
}


// PREDEFINED MODELS
// [[Rcpp::export]]
Rcpp::List run_algorithm(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg(x0, theta, S, tout, lower, upper, tau, rates_function));
}
std::tuple<arma::vec,bool,arma::vec> run_algorithm2(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg2(x0, theta, S, tout, lower, upper, tau, rates_function));
}
std::tuple<arma::vec,bool,arma::vec> run_algorithm_debug(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau, DebugRB &debug){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_debug(x0, theta, S, tout, lower, upper, tau, rates_function,debug));
  // return(gillespie_alg_debug2(x0, theta, S, tout, rates_function));
}
Rcpp::List run_algorithm_frac(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_frac(x0, theta, S, tout, rates_function));
}

std::tuple<arma::vec> run_algorithm_frac2(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_frac2(x0, theta, S, tout, rates_function));
}
std::tuple<arma::vec> run_algorithm_frac_debug(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout, DebugFrac &debug){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_frac_debug(x0, theta, S, tout, rates_function,debug));
}

// [[Rcpp::export]]
Rcpp::List sim_data(std::string const &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  return(gillespie_alg_entire(x0, theta, S, tout, lower, upper, tau, rates_function));
}

// [[Rcpp::export]]
Rcpp::List sim_data_frac(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout){
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
arma::mat cexpm(const arma::mat &mat){
  return(SS_exp_Q(Rcpp::NumericMatrix(Rcpp::wrap(mat)),10e-15,false));
  // return(arma::expmat(mat));
}
// [[Rcpp::export]]
arma::mat cexpm2(const arma::mat &mat){
  // return(SS_exp_Q(Rcpp::NumericMatrix(Rcpp::wrap(mat)),10e-15,false));
  return(arma::expmat(mat));
}

// [[Rcpp::export]]
double get_f(const arma::vec &xt, const arma::vec &yt){
  // Inexact Uniform:
  int box_halfwidth = 4;
  if(arma::sum(arma::abs(xt-yt))<=box_halfwidth){
    return(1/pow(2*box_halfwidth+1,xt.size()));
  }
  return(0);

  // Inexact Missing Prey:  
  // if(fabs(xt(1)-yt(1))<0.5){
  //   return(1);
  // }
  // return(0);

  // Inexact Gaussian Noise
  // arma::vec zt = xt-yt;
  // int d = zt.size();
  // double f = 1;
  // double noise_sd = 0.01;
  // for(int i=0;i<d;++i){
  //   f = f*R::dnorm(zt(i),0,noise_sd,0);
  // }
  // return(f);
  // Exact obs BDI:
  // if(fabs(xt(0)-yt(0))<0.5){
  //   return(1);
  // }
  // return(0);
  // Exact obs LV:
  // if(arma::sum(arma::abs(xt-yt))<=0.5){
  //   return(1);
  // }
  return(0);
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
  double est = hitObs*(1-inS) + P[state_to_index(xttau_data,lower,upper),state_to_index(xt_data,lower,upper)];
  return(est);
}

double get_estimate_noise_P(const std::tuple<arma::vec,bool,arma::vec> &estimate, const arma::vec &lower, const arma::vec &upper, const arma::mat &P, const arma::vec &obs, const arma::vec &fP){  
  arma::vec xttau_data = std::get<0>(estimate);
  bool inS = std::get<1>(estimate); 
  arma::vec xt_data = std::get<2>(estimate);
  int tau_index = state_to_index(xttau_data,lower,upper);
  double est = get_f(xt_data,obs)*(1-inS) + fP(tau_index);
  return(est);
}
//[[Rcpp::export]]
double get_estimate_noise_P2(const Rcpp::List &estimate, const arma::vec &lower, const arma::vec &upper, const arma::mat &P, const arma::vec &obs, const arma::vec &fP){  
  arma::vec xttau_data = estimate["xttau_data"];
  bool inS = estimate["inS"]; 
  arma::vec xt_data = estimate["xt_data"];
  int tau_index = state_to_index(xttau_data,lower,upper);
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

double get_estimate_frac_noise(const std::tuple<arma::vec> &estimate, const arma::vec &obs){  
  arma::vec xt_data = std::get<0>(estimate);
  return(get_f(xt_data,obs));
}

//[[Rcpp::export]]
double RB(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau,
                  const arma::vec &obs, const int &M){
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

uint64_t timeSinceEpochNanosec() {
  using namespace std::chrono;
  return duration_cast<nanoseconds>(system_clock::now().time_since_epoch()).count();
}

//[[Rcpp::export]]
std::string vec_to_string(const arma::vec &vec){
  std::ostringstream oss;
  oss.precision(5);
  for(int i=0; i<vec.n_elem;++i){
    oss<<vec[i]<<",";
  }
  return(oss.str());
}

void add_to_map(const std::string &key, const double &weight, std::unordered_map<std::string,double> &map){
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
arma::mat sample_x(const arma::mat &resample_vec_xt, const arma::vec &weights, const int &num_species, const int &M){
  int xt_size = weights.size(); 
  double normalise_constant = 0;
  for(int i=0; i<xt_size; ++i){
    normalise_constant += weights(i);
  }
  // std::random_device rd;
  // std::mt19937 gen(rd());
  // std::uniform_real_distribution<> dis(0.0,normalise_constant);
  arma::mat next_x(M,num_species);
  double rand_num;
  double cum_sum;
  for(int i=0;i<M;++i){
    // rand_num = dis(gen);
    rand_num = R::runif(0,normalise_constant);
    cum_sum = 0;
    for(int j=0;j<xt_size;++j){
      cum_sum += weights(j);
      if(rand_num<cum_sum){
        next_x.row(i) = resample_vec_xt.row(j);
        break;      
      }
    }
  }  
  return(next_x);
}

//[[Rcpp::export]]
arma::vec sample_x2(const arma::vec &x, const int &size, const bool &replace, const arma::vec &probs){
  return(Rcpp::RcppArmadillo::sample(x, size, replace, probs));
}

std::pair<arma::mat,arma::vec> sample_x_gibbs(std::vector<std::pair<std::pair<arma::vec,int>,double>> &resample_vec, const int &M, Swarm &swarm){
  int xt_size = resample_vec.size(); 
  double normalise_constant = 0;
  int num_species = resample_vec[0].first.first.size();

  for(const auto &entry : resample_vec){
    normalise_constant += entry.second;
  }
  std::vector<double> probabilities;
  for(auto &entry : resample_vec){
    entry.second /= normalise_constant;
    if(probabilities.empty()){
      probabilities.push_back(entry.second);
    }
    else{
      probabilities.push_back(entry.second+probabilities.back());
    }
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0,1.0);
  arma::mat next_x(M,num_species);
  std::vector<string> samples;
  samples.reserve(M);
  arma::vec parent_indices(M);
  std::pair<arma::vec,int> state;
  double probability;
  double rand_num;
  for(int i=0;i<M;++i){
    rand_num = dis(gen);
    auto it = std::lower_bound(probabilities.begin(),probabilities.end(),rand_num);
    int index = std::distance(probabilities.begin(),it); 
    next_x.row(i) = resample_vec[index].first.first.t();
    // parent_indices(i) = resample_vec[index].first.second;
    // parent_indices(i) = index;
    parent_indices(i) = i;
    state = std::make_pair(next_x.row(i).t(),parent_indices(i));
    probability = probabilities[index];
    swarm.add_state(state,probability);
  }  
  swarm.push_state();
  return(make_pair(next_x,parent_indices));
}

std::pair<arma::mat,double> RB_noise(const std::string &str, const arma::mat &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau,
                  const arma::vec &obs, const int &M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  arma::mat Q = get_coffin_matrix(str,lower,upper,theta);
  arma::mat P = cexpm(Q*tau);
  arma::vec f = get_ft(lower,upper,obs,P);
  int tau_index;
  std::tuple<arma::vec,bool,arma::vec> estimate;
  arma::mat resample_vec_xt(Q.n_rows-1+M,lower.size());
  arma::vec weights(Q.n_rows-1+M);
  double estimator = 0;
  arma::vec xttau_data(lower.size());
  arma::vec xt_data(lower.size());
  bool inS;
  arma::vec x;
  int resample_counter = Q.n_rows-1;
  int index = 0;
  for(int i=0;i<Q.n_rows-1;++i){
    resample_vec_xt.row(i) = index_to_state(i,lower,upper).t();
  }
  arma::mat Pf = P*f;
  for(int i=0;i<M;++i){
    x = x0.row(i).t();
    estimate = run_algorithm2(str,x,theta, tout, lower, upper, tau);
    estimator += get_estimate_noise_P(estimate,lower,upper,P,obs,Pf);
    xttau_data = std::get<0>(estimate);
    inS = std::get<1>(estimate);
    xt_data = std::get<2>(estimate);
    weights(resample_counter) = get_f(xt_data,obs)*(1-inS);
    resample_vec_xt.row(resample_counter) = xt_data.t();
    index = state_to_index(xttau_data,lower,upper);
    for(int j=0;j<Q.n_rows-1;++j){
      weights(j) += f(j)*P(index,j);
    }
    resample_counter += 1;
  }
  arma::vec s(weights.size());
  arma::mat next_x(M,obs.size());
  bool weights_check = false;
  for(int i=0;i<weights.size();++i){
    s(i)=i;
    if(weights(i)>0){
      weights_check = true;
    }
  }
  if(!weights_check){
    resample_counter = Q.n_rows-1;
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(resample_counter);
      resample_counter += 1;      
    }
  }
  else{
    arma::vec samples = sample_x2(s,M,true,weights);
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(samples(i));
    }
  }
  // arma::mat next_x = sample_x(resample_vec_xt,weights,obs.size(),M);
  std::pair<arma::mat,double> result = std::make_pair(next_x, estimator/M);
  return(result);
}

DebugRB RB_noise_debug(const std::string &str, const arma::mat &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau,
                  const arma::vec &obs, const int &M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  uint64_t t1 = timeSinceEpochNanosec();
  arma::mat Q = get_coffin_matrix(str,lower,upper,theta);
  arma::mat P = cexpm(Q*tau);
  uint64_t t2 = timeSinceEpochNanosec();
  arma::vec f = get_ft(lower,upper,obs,P);
  int tau_index;
  std::tuple<arma::vec,bool,arma::vec> estimate;
  arma::mat resample_vec_xt(Q.n_rows-1+M,lower.size());
  arma::vec weights(Q.n_rows-1+M);
  double estimator = 0;
  bool hitS;
  arma::vec xttau_data(lower.size());
  bool inS;
  arma::vec xt_data(lower.size());
  arma::vec inS_vec(M);
  arma::vec hitS_vec(M);
  arma::mat xt_data_mat(M,x0.n_cols);
  arma::mat xttau_data_mat(M,x0.n_cols);
  arma::vec x;
  DebugRB result;

  int resample_counter = Q.n_rows-1;
  int index = 0;
    for(int i=0;i<Q.n_rows-1;++i){
    resample_vec_xt.row(i) = index_to_state(i,lower,upper).t();
  }
  arma::mat Pf = P*f;
  arma::vec particle_weights(M);
  uint64_t before_gillespie;
  uint64_t after_gillespie;
  uint64_t after_weight;
  uint64_t gillespie_diff = 0.0;
  uint64_t weight_diff = 0.0;
  for(int i=0;i<M;++i){
    x = x0.row(i).t();
    before_gillespie = timeSinceEpochNanosec();
    estimate = run_algorithm_debug(str,x,theta, tout, lower, upper, tau,result);
    after_gillespie = timeSinceEpochNanosec();
    particle_weights(i) = get_estimate_noise_P(estimate,lower,upper,P,obs,Pf);
    after_weight = timeSinceEpochNanosec();
    gillespie_diff += after_gillespie-before_gillespie;
    weight_diff += after_weight-after_gillespie;

    estimator += particle_weights(i);

    xttau_data = std::get<0>(estimate);
    inS = std::get<1>(estimate);
    xt_data = std::get<2>(estimate);
    weights(resample_counter) = get_f(xt_data,obs)*(1-inS);
    resample_vec_xt.row(resample_counter) = xt_data.t();
    index = state_to_index(xttau_data,lower,upper);
    for(int j=0;j<Q.n_rows-1;++j){
      weights(j) += f(j)*P(index,j);
    }
    resample_counter += 1;
    hitS = true;
    for(int j=0;j<xttau_data.size();++j){
      if((xttau_data[j] > upper[j]) || (lower[j] > xttau_data[j])){
        hitS = false;
      }
    }
    inS_vec(i) = inS;
    hitS_vec(i) = hitS;
    xt_data_mat.row(i) = xt_data.t();
    xttau_data_mat.row(i) = xttau_data.t();
  }
  uint64_t before_resample = timeSinceEpochNanosec();
  arma::vec s(weights.size());
  arma::mat next_x(M,obs.size());
  bool weights_check = false;
  for(int i=0;i<weights.size();++i){
    s(i)=i;
    if(weights(i)>0){
      weights_check = true;
    }
  }
  if(!weights_check){
    resample_counter = Q.n_rows-1;
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(resample_counter);
      resample_counter += 1;      
    }
  }
  else{
    arma::vec samples = sample_x2(s,M,true,weights);
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(samples(i));
    }
  }
  // arma::mat next_x = sample_x(resample_vec_xt,weights,obs.size(),M);
  uint64_t t3 = timeSinceEpochNanosec();
  uint64_t exp_time = t2-t1;
  uint64_t resample_time = t3-before_resample;
  uint64_t entire_time = t3-t1;
  result.x=next_x;
  result.log_prob=estimator/M;
  result.inS=inS_vec;
  result.hitS=hitS_vec;
  result.xt_data=xt_data_mat;
  result.xttau_data=xttau_data_mat;
  result.presample=resample_vec_xt;
  result.weights=weights;
  result.particle_weights=particle_weights;
  result.exp_time=exp_time;
  result.entire_time=entire_time;
  result.gillespie_time=gillespie_diff;
  result.weight_time=weight_diff;
  result.resample_time=resample_time;
  return(result);
}

std::pair<std::pair<arma::mat,arma::vec>,double> RB_noise_gibbs(const std::string &str, const std::pair<arma::mat,arma::vec> &x0_pair, const arma::vec &theta, const double &tout,
                  const arma::vec &lower, const arma::vec &upper, const double &tau,
                  const arma::vec &obs, const int &M, 
                  Swarm &swarm, arma::vec &numInS, const int &obs_index){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  arma::mat Q = get_coffin_matrix(str,lower,upper,theta).t();
  arma::mat P = cexpm(Q*tau);
  arma::vec f = get_ft(lower,upper,obs,P);
  int tau_index;
  std::tuple<arma::vec,bool,arma::vec> estimate;
  std::vector<std::pair<std::pair<arma::vec,int>,double>> resample_vec;
  double estimator = 0;
  arma::vec xttau_data;
  bool inS;
  arma::vec xt_data;
  arma::vec x;
  int parent_index;
  arma::vec fP;
  for(int i=0;i<M;++i){
    x = x0_pair.first.row(i).t();
    parent_index = x0_pair.second(i);
    estimate = run_algorithm2(str,x,theta, tout, lower, upper, tau);
    estimator += get_estimate_noise_P(estimate,lower,upper,P,obs,(f.t()*P).t());
    inS = std::get<1>(estimate);
    xt_data = std::get<2>(estimate);
    numInS[obs_index] = numInS[obs_index] + inS;
    resample_vec.push_back(make_pair(std::make_pair(xt_data,parent_index),get_f(xt_data,obs)*(1-inS)));
    fP = f%P.col(state_to_index(xttau_data,lower,upper));
    for(int j=0;j<fP.size()-1;++j){
      resample_vec.push_back(make_pair(make_pair(index_to_state(j,lower,upper),parent_index),fP(j)));
    }
  }
  std::pair<arma::mat,arma::vec> next_x = sample_x_gibbs(resample_vec,M,swarm);
  std::pair<std::pair<arma::mat,arma::vec>,double> result = std::make_pair(next_x, estimator/M);
  return(result);
}

//[[Rcpp::export]]
double RB_list(const std::string &str, const arma::vec &x0, const arma::vec &theta, const arma::vec &tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, const arma::vec &tau_list,
                  const arma::mat &obs_list, const int &M){
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
    lower = lower_list.row(i).t();
    upper = upper_list.row(i).t();
    obs = obs_list.row(i).t();
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
    tau = tau_list(i);
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
Rcpp::List RB_particle_filter_debug(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  Rcpp::List debug_list(tout_list.size());
  arma::vec lower(num_species);
  arma::vec upper(num_species);
  arma::vec obs(num_species);
  double tout;
  double tau;
  arma::mat x = arma::repmat(x0.t(),M,1);
  double prev_time = 0;
  double log_prob = 0;
  // std::tuple<arma::mat,double,arma::vec,arma::vec,arma::mat,arma::mat,arma::mat,arma::vec,arma::vec,uint64_t,uint64_t> result;
  DebugRB result;
  for(int i=0; i<num_obs; ++i){
    lower = lower_list.row(i).t();
    upper = upper_list.row(i).t();
    obs = obs_list.row(i).t();
    tout = tout_list(i)-prev_time;
    tau = tau_list(i);
    result = RB_noise_debug(str,x,theta,tout,lower,upper,tau,obs,M);
    log_prob += log(result.log_prob);
    x = result.x; // copy
    prev_time = tout_list(i);
    Rcpp::List tupleList = result.get_tuple();
    debug_list[i] = tupleList;
    if(isinf(log_prob)){
      return(debug_list);
    }
  }
  return(debug_list);
}

//[[Rcpp::export]]
double frac(const std::string &str, const arma::vec &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &obs, const int &M){
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

std::pair<arma::mat,double> frac_noise(const std::string &str, const arma::mat &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &obs, const int &M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  std::tuple<arma::mat> estimate;
  arma::mat resample_vec_xt(M,obs.size());
  arma::vec weights(M);
  double estimator = 0;
  arma::vec x;
  arma::vec xt_data(obs.size());
  for(int i=0;i<M;++i){
    x = x0.row(i).t();
    estimate = run_algorithm_frac2(str,x,theta, tout);
    xt_data = std::get<0>(estimate);
    estimator += get_estimate_frac_noise(estimate,obs);
    resample_vec_xt.row(i) = xt_data.t();
    weights(i) = get_f(xt_data,obs);
  }
    arma::vec s(weights.size());
  arma::mat next_x(M,obs.size());
  bool weights_check = false;
  for(int i=0;i<weights.size();++i){
    s(i)=i;
    if(weights(i)>0){
      weights_check = true;
    }
  }
  if(!weights_check){
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(i);
    }
  }
  else{
    arma::vec samples = sample_x2(s,M,true,weights);
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(samples(i));
    }
  }
  // arma::mat next_x = sample_x(resample_vec_xt,weights,obs.size(),M);
  std::pair<arma::mat,double> result = std::make_pair(next_x, estimator/M);
  return(result);
}

DebugFrac frac_noise_debug(const std::string &str, const arma::mat &x0, const arma::vec &theta, const double &tout,
                  const arma::vec &obs, const int &M){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  arma::mat S = get_S(str);
  std::tuple<arma::mat> estimate;
  arma::mat resample_vec_xt(M,obs.size());
  arma::vec weights(M);
  double estimator = 0;
  arma::vec x;
  arma::vec xt_data(obs.size());
  uint64_t t1 = timeSinceEpochNanosec(); 
  uint64_t before_gillespie; 
  uint64_t after_gillespie; 
  uint64_t after_weight; 
  uint64_t gillespie_diff = 0.0; 
  uint64_t weight_diff = 0.0; 

  for(int i=0;i<M;++i){
    x = x0.row(i).t();
    before_gillespie = timeSinceEpochNanosec();    
    estimate = run_algorithm_frac2(str,x,theta, tout);
    after_gillespie = timeSinceEpochNanosec();
    xt_data = std::get<0>(estimate);
    estimator += get_estimate_frac_noise(estimate,obs);
    after_weight = timeSinceEpochNanosec();
    gillespie_diff += after_gillespie-before_gillespie;
    weight_diff += after_weight-after_gillespie;

    resample_vec_xt.row(i) = xt_data.t();
    weights(i) = get_f(xt_data,obs);
  }
  uint64_t before_resample = timeSinceEpochNanosec();
  arma::vec s(weights.size());
  arma::mat next_x(M,obs.size());
  bool weights_check = false;
  for(int i=0;i<weights.size();++i){
    s(i)=i;
    if(weights(i)>0){
      weights_check = true;
    }
  }
  if(!weights_check){
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(i);
    }
  }
  else{
    arma::vec samples = sample_x2(s,M,true,weights);
    for(int i=0;i<M;++i){
      next_x.row(i) = resample_vec_xt.row(samples(i));
    }
  }
  // arma::mat next_x = sample_x(resample_vec_xt,weights,obs.size(),M);
  uint64_t t2 = timeSinceEpochNanosec(); 
  uint64_t resample_time = t2-before_resample; 
  uint64_t entire_time = t2-t1; 
  DebugFrac result;
  result.x = next_x;
  result.log_prob = estimator/M;
  result.presample = resample_vec_xt;
  result.weights = weights;
  result.entire_time = entire_time;
  result.gillespie_time=gillespie_diff;
  result.weight_time=weight_diff;
  result.resample_time=resample_time;
  return(result);
}

std::pair<std::pair<arma::mat,arma::vec>,double> frac_noise_gibbs(const std::string &str, const std::pair<arma::mat,arma::vec> &x0_pair, const arma::vec &theta, const double &tout,
                  const arma::vec &obs, const int &M, 
                  Swarm &swarm){
  std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function = get_rate_function(str);
  std::tuple<arma::vec> estimate;
  std::vector<std::pair<std::pair<arma::vec,int>,double>> resample_vec;
  double estimator = 0;
  arma::vec xt_data;
  arma::vec x;
  int parent_index;
  for(int i=0;i<M;++i){
    x = x0_pair.first.row(i).t();
    parent_index = x0_pair.second(i);
    xt_data = std::get<0>(estimate);
    estimate = run_algorithm_frac2(str,x,theta, tout);
    estimator += get_estimate_frac_noise(estimate,obs);
    resample_vec.push_back(make_pair(std::make_pair(xt_data,parent_index),get_f(xt_data,obs)));
  }
  std::pair<arma::mat,arma::vec> next_x = sample_x_gibbs(resample_vec,M,swarm);
  std::pair<std::pair<arma::mat,arma::vec>,double> result = std::make_pair(next_x, estimator/M);
  return(result);
}


//[[Rcpp::export]]
double frac_list(const std::string &str, const arma::vec &x0, const arma::vec &theta, const arma::vec &tout_list,
                  const arma::mat &obs_list, const int &M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec obs(num_species);
  double tout;
  arma::vec x = x0; // copy
  double prev_time = 0;
  double log_prob = 0;
  for(int i=0; i<num_obs; ++i){
    obs = obs_list.row(i).t();
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
  arma::mat x = arma::repmat(x0.t(),M,1);
  double prev_time = 0;
  double log_prob = 0;
  std::pair<arma::mat,double> result;
  for(int i=0; i<num_obs; ++i){
    obs = obs_list.row(i).t();
    tout = tout_list(i)-prev_time;
    result = frac_noise(str,x,theta,tout,obs,M);
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
Rcpp::List frac_particle_filter_debug(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &obs_list, int M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec obs(num_species);
  double tout;
  arma::mat x = arma::repmat(x0.t(),M,1);
  double prev_time = 0;
  double log_prob = 0;
  Rcpp::List debug_list(tout_list.size());
  DebugFrac result;
  for(int i=0; i<num_obs; ++i){
    obs = obs_list.row(i).t();
    tout = tout_list(i)-prev_time;
    result = frac_noise_debug(str,x,theta,tout,obs,M);
    Rcpp::List tupleList = result.get_tuple();
    debug_list[i] = tupleList;
    log_prob += log(result.log_prob);
    if(isinf(log_prob)){
      return(debug_list);
    }
    x = result.x; // copy
    prev_time = tout_list(i);
  }
  return(debug_list);
}

//[[Rcpp::export]]
Rcpp::List RB_particle_filter_gibbs(std::string str, const arma::vec &x0, const arma::vec &theta, arma::vec tout_list,
                  const arma::mat &lower_list, const arma::mat &upper_list, arma::vec tau_list,
                  const arma::mat &obs_list, int M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec lower(num_species);
  arma::vec upper(num_species);
  arma::vec obs(num_species);
  double tout;
  double tau;
  arma::vec parent_index = arma::vec(M);
  std::pair<arma::mat,arma::vec> x = std::make_pair(arma::repmat(x0.t(),M,1),parent_index);
  double prev_time = 0;
  double log_prob = 0;
  std::pair<std::pair<arma::mat,arma::vec>,double> result;
  Swarm swarm;
  arma::vec numInS(num_obs);
  for(int i=0; i<num_obs; ++i){
    lower = lower_list.row(i).t();
    upper = upper_list.row(i).t();
    obs = obs_list.row(i).t();
    tout = tout_list(i)-prev_time;
    tau = tau_list[i];
    result = RB_noise_gibbs(str,x,theta,tout,lower,upper,tau,obs,M,swarm,numInS,i);
    log_prob += log(result.second);
    if(isinf(log_prob)){
      return(log_prob);
    }
    x = result.first; // copy
    prev_time = tout_list(i);
  }
  arma::mat path = swarm.get_path();
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("path")=path,Rcpp::Named("ll")=log_prob, Rcpp::Named("numInS")=numInS);
  return(results);
}

//[[Rcpp::export]]
Rcpp::List frac_particle_filter_gibbs(const std::string &str, const arma::vec &x0, const arma::vec &theta, const arma::vec &tout_list,
                  const arma::mat &obs_list, const int &M){
  int num_obs = tout_list.n_elem;
  int num_species = x0.n_elem;
  arma::vec obs(num_species);
  double tout;
  arma::vec parent_index = arma::vec(M);
  std::pair<arma::mat,arma::vec> x = std::make_pair(arma::repmat(x0.t(),M,1),parent_index);
  double prev_time = 0;
  double log_prob = 0;
  std::pair<std::pair<arma::mat,arma::vec>,double> result;
  Swarm swarm;
  for(int i=0; i<num_obs; ++i){
    obs = obs_list.row(i).t();
    tout = tout_list(i)-prev_time;
    result = frac_noise_gibbs(str,x,theta,tout,obs,M,swarm);
    log_prob += log(result.second);
    if(isinf(log_prob)){
      return(log_prob);
    }
    x = result.first; // copy
    prev_time = tout_list(i);
  }
  arma::mat path = swarm.get_path();
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("path")=path,Rcpp::Named("ll")=log_prob);
  return(results);
}

//[[Rcpp::export]]
arma::vec test_RB(const std::string &str,const arma::vec &x,const arma::vec &theta, const double &tout,const arma::vec &lower, const arma::vec &upper, const double &tau, const int &M){
  u_int64_t before_gillespie;
  u_int64_t after_gillespie;
  u_int64_t total_gillespie = 0.0;
  arma::mat Q = get_coffin_matrix(str,lower,upper,theta);
  arma::mat P = cexpm(Q*tau);
  std::tuple<arma::vec,bool,arma::vec> estimate;
  DebugRB debug;
  arma::vec results(M);
  for(int i=0;i<M;++i){
    before_gillespie = timeSinceEpochNanosec();    
    estimate = run_algorithm_debug(str,x,theta, tout, lower,upper,tau,debug);
    after_gillespie = timeSinceEpochNanosec();
    total_gillespie = after_gillespie-before_gillespie;
    results(i) = total_gillespie;
  }
  return(results);
}
//[[Rcpp::export]]
arma::vec test_frac(const std::string &str,const arma::vec &x,const arma::vec &theta, const double &tout,const arma::vec &lower, const arma::vec &upper, const double &tau, const int &M){
  u_int64_t before_gillespie;
  u_int64_t after_gillespie;
  u_int64_t total_gillespie = 0.0;
  // Rcpp::NumericMatrix Q = Rcpp::NumericMatrix(Rcpp::wrap(get_coffin_matrix(str,lower,upper,theta)));
  // P = cexpm(Q*tau);
  // P = SS_exp_Q(Q*tau,10e-10,false);    
  std::tuple<arma::mat> estimate;
  DebugFrac debug;
  arma::vec results(M);
  for(int i=0;i<M;++i){
    before_gillespie = timeSinceEpochNanosec();    
    estimate = run_algorithm_frac_debug(str,x,theta, tout,debug);
    after_gillespie = timeSinceEpochNanosec();
    total_gillespie = after_gillespie-before_gillespie;
    results(i)=total_gillespie;
  }
  return(results);
}