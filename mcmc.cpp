#include "mcmc.h"

//[[Rcpp::export]]
Rcpp::List mcmc_RB(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const int &iterations, const int &num_particles){
  int num_params = thetas.size();
  Rcpp::NumericMatrix posterior(iterations,num_params);
  Rcpp::NumericVector ll(iterations);
  Rcpp::NumericVector log_theta = log(thetas);
  Rcpp::NumericVector log_psi;
  double psi_ll;
  double psi_posterior;
  double alpha;
  Rcpp::NumericVector unif = Rcpp::runif(iterations,0.0,1.0);
  // arma::mat sigma = arma::mat("0.4400080, 0.1306354, -0.2482663; 0.1306354, 0.1551602,  0.1469131; -0.2482663, 0.1469131,  0.7858287");
  arma::mat sigma = arma::mat("0.4738750, 0.1448164, -0.1204549;0.1448164, 0.1705808,  0.1844356; -0.1204549, 0.1844356,  0.5225753");
  Rcpp::NumericMatrix offsets = Rcpp::wrap(mvrnormArma(iterations,sigma));
  double theta_ll = get_likelihood_RB(model_name,start,exp(log_theta),obs_list,tout_list,num_particles);
  double theta_posterior = get_prior(log_theta)+theta_ll;
  double acceptance=0;
  for(int i=0;i<iterations;++i){
    // log_psi = propose_psi(log_theta); 
    log_psi = log_theta + offsets(i,Rcpp::_)*1.374;
    std::cout<<"Iteration:"<<i<<", Acceptance:"<<acceptance/(i+1)<<", Proposal"<<log_psi<<std::endl;;
    psi_ll = get_likelihood_RB(model_name,start,exp(log_psi),obs_list,tout_list,num_particles);
    psi_posterior = get_prior(log_psi)+psi_ll;
    // alpha = psi_ll-theta_ll;
    alpha = psi_posterior-theta_posterior;    
    if(isnan(alpha)){
      alpha = -DBL_MAX;
    }
    // std::cout<<"Iteration:"<<i<<", Theta:"<<(log_theta)<<", ThetaLL:"<<theta_ll<<" | "<<"Psi:"<<(log_psi)<<", PsiLL:"<<psi_ll<<", Alpha: "<<alpha<<std::endl;;
    if(log(unif[i])<alpha){
      // std::cout<<", Accepted."<<std::endl;
      log_theta = clone(log_psi);
      theta_ll = psi_ll;
      theta_posterior = psi_posterior;
      acceptance+=1;
    }
    else{
      // std::cout<<", Rejected."<<std::endl;
    }
    posterior(i,Rcpp::_) = log_theta;
    ll(i) = theta_ll;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("log_posterior")=posterior,Rcpp::Named("ll")=ll,Rcpp::Named("acceptance_rate")=acceptance/iterations);
  return(results);
}
//[[Rcpp::export]]
Rcpp::List mcmc_frac(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const int &iterations, const int &num_particles){
  int num_params = thetas.size();
  Rcpp::NumericMatrix posterior(iterations,num_params);
  Rcpp::NumericVector ll(iterations);
  Rcpp::NumericVector log_theta = log(thetas);
  Rcpp::NumericVector log_psi;
  double psi_ll;
  double psi_posterior;
  double alpha;
  Rcpp::NumericVector unif = Rcpp::runif(iterations,0.0,1.0);
  // arma::mat sigma = arma::mat("1,0,0; 0,1,0; 0,0,1");
  // arma::mat sigma = arma::mat("0.4400080, 0.1306354, -0.2482663; 0.1306354, 0.1551602,  0.1469131; -0.2482663, 0.1469131,  0.7858287");
  arma::mat sigma = arma::mat("0.4738750, 0.1448164, -0.1204549; 0.1448164, 0.1705808,  0.1844356; -0.1204549, 0.1844356,  0.5225753");
  Rcpp::NumericMatrix offsets = Rcpp::wrap(mvrnormArma(iterations,sigma));
  double theta_ll = get_likelihood_frac(model_name,start,exp(log_theta),obs_list,tout_list,num_particles);
  double theta_posterior = get_prior(log_theta)+theta_ll;
  double acceptance=0;
  for(int i=0;i<iterations;++i){
    // log_psi = propose_psi(log_theta);    
    // log_psi = log_theta + offsets(i,Rcpp::_)*1.5;
    log_psi = log_theta + offsets(i,Rcpp::_)*1.374; 
    std::cout<<"Iteration:"<<i<<", Acceptance:"<<acceptance/(i+1)<<", Proposal"<<log_psi<<std::endl;;
    psi_ll = get_likelihood_frac(model_name,start,exp(log_psi),obs_list,tout_list,num_particles);
    psi_posterior = get_prior(log_psi)+psi_ll;
    // alpha = psi_ll-theta_ll;
    alpha = psi_posterior-theta_posterior;    
    if(isnan(alpha)){
      alpha = -DBL_MAX;
    }
    // std::cout<<"Iteration:"<<i<<", Theta:"<<(log_theta)<<", ThetaLL:"<<theta_ll<<" | "<<"Psi:"<<(log_psi)<<", PsiLL:"<<psi_ll<<", Alpha: "<<alpha<<std::endl;;
    if(log(unif[i])<alpha){
      // std::cout<<", Accepted."<<std::endl;
      log_theta = clone(log_psi);
      theta_ll = psi_ll;
      theta_posterior = psi_posterior;
      acceptance+=1;
    }
    else{
      // std::cout<<", Rejected."<<std::endl;
    }
    posterior(i,Rcpp::_) = log_theta;
    ll(i) = theta_ll;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("log_posterior")=posterior,Rcpp::Named("ll")=ll,Rcpp::Named("acceptance_rate")=acceptance/iterations);
  return(results);
}

//[[Rcpp::export]]
Rcpp::NumericVector propose_psi(const Rcpp::NumericVector &theta){
  Rcpp::NumericVector offset = Rcpp::rnorm(3,0.0,0.1);
  Rcpp::NumericVector psi = theta+offset;
  return(psi);
}

//[[Rcpp::export]]
double get_likelihood_RB(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &theta, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const int &num_particles){
  double goal = sqrt(0.95);
  Rcpp::NumericVector tau_list = get_tau_list(model_name,start,theta,obs_list,tout_list,goal);
  Rcpp::List box = get_box_brownian_fast_list(model_name,theta,tout_list,tau_list,start,obs_list,goal);
  Rcpp::NumericMatrix lower_list = box["lower"];
  Rcpp::NumericMatrix upper_list = box["upper"];
  double likelihood = RB_list(model_name, start, theta, tout_list, lower_list, upper_list, tau_list, obs_list, num_particles);
  return(likelihood);
}

//[[Rcpp::export]]
double get_likelihood_frac(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &theta, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const int &num_particles){
  double likelihood = frac_list(model_name, start, theta, tout_list, obs_list, num_particles);
  return(likelihood);
}

// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(n, ncols);
  //  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
   return(Y * arma::chol(sigma));
}

//[[Rcpp::export]]
double get_prior(const Rcpp::NumericVector log_theta){
  int num_params = log_theta.size();
  // Rcpp::NumericVector mus = {0.2,0.3,10};
  Rcpp::NumericVector mus = {-1.6094,-1.204,2.303};
  Rcpp::NumericVector sigmas = {1,1,1};
  double log_prior = 0;
  for(int i=0;i<num_params;++i){
    log_prior += -((log_theta[i]-mus[i])*(log_theta[i]-mus[i]))/(2*sigmas[i]);
  }
  return(log_prior);
}