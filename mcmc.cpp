#include "mcmc.h"

//[[Rcpp::export]]
Rcpp::List mcmc_RB(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::mat &obs_list, const arma::mat &lower_list, const arma::mat &upper_list, 
const arma::vec &tout_list, const int &iterations, const int &num_particles, const arma::mat &sigma){
  int num_params = thetas.size();
  arma::mat posterior(iterations,num_params);
  arma::vec ll(iterations);
  arma::vec transformed_theta = transform(thetas);
  arma::vec transformed_psi;
  double psi_ll;
  double psi_posterior;
  double alpha;
  arma::vec unif = Rcpp::runif(iterations,0.0,1.0);
  arma::mat offsets = mvrnormArma(iterations,sigma);
  double theta_ll = get_likelihood_RB(model_name,start,inv_transform(transformed_theta),obs_list,lower_list,upper_list,tout_list,num_particles);
  double theta_posterior = get_prior(transformed_theta)+theta_ll;
  double acceptance=0;
  double scale = sqrt(2.38*2.38/num_params);
  for(int i=0;i<iterations;++i){
    // transformed_psi = propose_psi(transformed_theta); 
    transformed_psi = transformed_theta + offsets.col(i)*scale;
    std::cout<<"Iteration:"<<i<<", Acceptance:"<<acceptance/(i+1)<<", Proposal"<<transformed_psi.t()<<std::endl;;
    psi_ll = get_likelihood_RB(model_name,start,inv_transform(transformed_psi),obs_list,lower_list,upper_list,tout_list,num_particles);
    psi_posterior = get_prior(transformed_psi)+psi_ll;
    // alpha = psi_ll-theta_ll;
    alpha = psi_posterior-theta_posterior;    
    if(isnan(alpha)){
      alpha = -DBL_MAX;
    }
    // std::cout<<"Iteration:"<<i<<", Theta:"<<(transformed_theta)<<", ThetaLL:"<<theta_ll<<" | "<<"Psi:"<<(transformed_psi)<<", PsiLL:"<<psi_ll<<", Alpha: "<<alpha<<std::endl;;
    if(log(unif[i])<alpha){
      // std::cout<<", Accepted."<<std::endl;
      transformed_theta = transformed_psi; // copy
      theta_ll = psi_ll;
      theta_posterior = psi_posterior;
      acceptance+=1;
    }
    else{
      // std::cout<<", Rejected."<<std::endl;
    }
    posterior.row(i) = transformed_theta.t();
    ll(i) = theta_ll;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("log_posterior")=posterior,Rcpp::Named("ll")=ll,Rcpp::Named("acceptance_rate")=acceptance/iterations);
  return(results);
}
//[[Rcpp::export]]
Rcpp::List mcmc_RB_marginal(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::mat &obs_list, const arma::mat &lower_list, const arma::mat &upper_list, 
const arma::vec &tout_list, const int &iterations, const int &num_particles){
  int num_params = thetas.size();
  arma::mat posterior(iterations,num_params);
  arma::vec ll(iterations);
  double psi_ll;
  double psi_posterior;
  double alpha;
  arma::vec unif = Rcpp::runif(iterations,0.0,1.0);
  Rcpp::List sim = get_likelihood_RB_gibbs(model_name,start,thetas,obs_list,lower_list,upper_list,tout_list,num_particles);
  double theta_ll = sim["ll"];
  double theta_posterior = theta_ll;
  double acceptance=0;
  for(int i=0;i<iterations;++i){
    std::cout<<"Iteration:"<<i<<", Acceptance:"<<acceptance/(i+1)<<std::endl;;
    sim = get_likelihood_RB_gibbs(model_name,start,thetas,obs_list,lower_list,upper_list,tout_list,num_particles);
    psi_ll = sim["ll"]; 
    psi_posterior = psi_ll;
    alpha = psi_posterior-theta_posterior;    
    if(isnan(alpha)){
      alpha = -DBL_MAX;
    }
    // std::cout<<"Iteration:"<<i<<", Theta:"<<(transformed_theta)<<", ThetaLL:"<<theta_ll<<" | "<<"Psi:"<<(transformed_psi)<<", PsiLL:"<<psi_ll<<", Alpha: "<<alpha<<std::endl;;
    if(log(unif[i])<alpha){
      // std::cout<<", Accepted."<<std::endl;
      theta_ll = psi_ll;
      theta_posterior = psi_posterior;
      acceptance+=1;
    }
    else{
      // std::cout<<", Rejected."<<std::endl;
    }
    posterior.row(i) = thetas.t();
    ll(i) = theta_ll;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("log_posterior")=posterior,Rcpp::Named("ll")=ll,Rcpp::Named("acceptance_rate")=acceptance/iterations);
  return(results);
}
//[[Rcpp::export]]
Rcpp::List mcmc_frac(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &iterations, const int &num_particles){
  int num_params = thetas.size();
  arma::mat posterior(iterations,num_params);
  arma::vec ll(iterations);
  arma::vec transformed_theta = transform(thetas);
  arma::vec transformed_psi;
  double psi_ll;
  double psi_posterior;
  double alpha;
  arma::vec unif = Rcpp::runif(iterations,0.0,1.0);
  // arma::mat sigma = arma::mat("1,0,0; 0,1,0; 0,0,1");
  // arma::mat sigma = arma::mat("0.4400080, 0.1306354, -0.2482663; 0.1306354, 0.1551602,  0.1469131; -0.2482663, 0.1469131,  0.7858287");
  arma::mat sigma = arma::mat("0.4738750, 0.1448164, -0.1204549; 0.1448164, 0.1705808,  0.1844356; -0.1204549, 0.1844356,  0.5225753");
  arma::mat offsets = mvrnormArma(iterations,sigma);
  double theta_ll = get_likelihood_frac(model_name,start,inv_transform(transformed_theta),obs_list,tout_list,num_particles);
  double theta_posterior = get_prior(transformed_theta)+theta_ll;
  double acceptance=0;
  double scale = sqrt(2.38*2.38/num_params);
  for(int i=0;i<iterations;++i){
    // transformed_psi = propose_psi(transformed_theta);    
    // transformed_psi = transformed_theta + offsets.row(i)*1.5;
    transformed_psi = transformed_theta + offsets.col(i)*scale; 
    std::cout<<"Iteration:"<<i<<", Acceptance:"<<acceptance/(i+1)<<", Proposal"<<transformed_psi<<std::endl;;
    psi_ll = get_likelihood_frac(model_name,start,inv_transform(transformed_psi),obs_list,tout_list,num_particles);
    psi_posterior = get_prior(transformed_psi)+psi_ll;
    // alpha = psi_ll-theta_ll;
    alpha = psi_posterior-theta_posterior;    
    if(isnan(alpha)){
      alpha = -DBL_MAX;
    }
    // std::cout<<"Iteration:"<<i<<", Theta:"<<(transformed_theta)<<", ThetaLL:"<<theta_ll<<" | "<<"Psi:"<<(transformed_psi)<<", PsiLL:"<<psi_ll<<", Alpha: "<<alpha<<std::endl;;
    if(log(unif[i])<alpha){
      // std::cout<<", Accepted."<<std::endl;
      transformed_theta = transformed_psi; // copy
      theta_ll = psi_ll;
      theta_posterior = psi_posterior;
      acceptance+=1;
    }
    else{
      // std::cout<<", Rejected."<<std::endl;
    }
    posterior.row(i) = transformed_theta;
    ll(i) = theta_ll;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("log_posterior")=posterior,Rcpp::Named("ll")=ll,Rcpp::Named("acceptance_rate")=acceptance/iterations);
  return(results);
}
//[[Rcpp::export]]
Rcpp::List mcmc_frac_marginal(const std::string &model_name, const arma::vec &start, const arma::vec &thetas, 
const arma::mat &obs_list, 
const arma::vec &tout_list, const int &iterations, const int &num_particles){
  int num_params = thetas.size();
  arma::mat posterior(iterations,num_params);
  arma::vec ll(iterations);
  double psi_ll;
  double psi_posterior;
  double alpha;
  arma::vec unif = Rcpp::runif(iterations,0.0,1.0);
  Rcpp::List sim = get_likelihood_frac_gibbs(model_name,start,thetas,obs_list,tout_list,num_particles);
  double theta_ll = sim["ll"];
  double theta_posterior = theta_ll;
  double acceptance=0;
  for(int i=0;i<iterations;++i){
    std::cout<<"Iteration:"<<i<<", Acceptance:"<<acceptance/(i+1)<<std::endl;;
    sim = get_likelihood_frac_gibbs(model_name,start,thetas,obs_list,tout_list,num_particles);
    psi_ll = sim["ll"];
    psi_posterior = psi_ll;
    alpha = psi_posterior-theta_posterior;    
    if(isnan(alpha)){
      alpha = -DBL_MAX;
    }
    // std::cout<<"Iteration:"<<i<<", Theta:"<<(transformed_theta)<<", ThetaLL:"<<theta_ll<<" | "<<"Psi:"<<(transformed_psi)<<", PsiLL:"<<psi_ll<<", Alpha: "<<alpha<<std::endl;;
    if(log(unif[i])<alpha){
      // std::cout<<", Accepted."<<std::endl;
      theta_ll = psi_ll;
      theta_posterior = psi_posterior;
      acceptance+=1;
    }
    else{
      // std::cout<<", Rejected."<<std::endl;
    }
    posterior.row(i) = thetas.t();
    ll(i) = theta_ll;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("log_posterior")=posterior,Rcpp::Named("ll")=ll,Rcpp::Named("acceptance_rate")=acceptance/iterations);
  return(results);
}

//[[Rcpp::export]]
arma::vec propose_psi(const arma::vec &theta){
  arma::vec offset = Rcpp::rnorm(3,0.0,0.1);
  arma::vec psi = theta+offset;
  return(psi);
}

//[[Rcpp::export]]
double get_likelihood_RB(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::mat &lower_list, const arma::mat &upper_list, const arma::vec &tout_list, const int &num_particles){
  double goal = 100;
  arma::vec tau_list = get_tau_list(model_name,start,theta,obs_list,tout_list,goal);
  // Rcpp::List box = get_box_brownian_fast_list(model_name,theta,tout_list,tau_list,start,obs_list,goal);
  // arma::mat lower_list = box["lower"];
  // arma::mat upper_list = box["upper"];
  // double likelihood = RB_list(model_name, start, theta, tout_list, lower_list, upper_list, tau_list, obs_list, num_particles);
  double likelihood = RB_particle_filter(model_name, start, theta, tout_list, lower_list, upper_list, tau_list, obs_list, num_particles);
  return(likelihood);
}

//[[Rcpp::export]]
double get_likelihood_frac(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &num_particles){
  // double likelihood = frac_list(model_name, start, theta, tout_list, obs_list, num_particles);
  double likelihood = frac_particle_filter(model_name, start, theta, tout_list, obs_list, num_particles);
  return(likelihood);
}

//[[Rcpp::export]]
Rcpp::List get_likelihood_RB_gibbs(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::mat &lower_list, const arma::mat &upper_list, const arma::vec &tout_list, const int &num_particles){
  double goal = 100;
  arma::vec tau_list = get_tau_list(model_name,start,theta,obs_list,tout_list,goal);
  // Rcpp::List box = get_box_brownian_fast_list(model_name,theta,tout_list,tau_list,start,obs_list,goal);
  // arma::mat lower_list = box["lower"];
  // arma::mat upper_list = box["upper"];
  // double likelihood = RB_list(model_name, start, theta, tout_list, lower_list, upper_list, tau_list, obs_list, num_particles);
  Rcpp::List result = RB_particle_filter_gibbs(model_name, start, theta, tout_list, lower_list, upper_list, tau_list, obs_list, num_particles);
  return(result);
}

//[[Rcpp::export]]
Rcpp::List get_likelihood_frac_gibbs(const std::string &model_name, const arma::vec &start, const arma::vec &theta, 
const arma::mat &obs_list, const arma::vec &tout_list, const int &num_particles){
  // double likelihood = frac_list(model_name, start, theta, tout_list, obs_list, num_particles);
  Rcpp::List result = frac_particle_filter(model_name, start, theta, tout_list, obs_list, num_particles);
  return(result);
}

// // [[Rcpp::export]]
// arma::mat mvrnormArma(int n, arma::mat sigma) {
//    int ncols = sigma.n_cols;
//    arma::mat Y = arma::randn(n, ncols);
//   //  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
//    return(Y * arma::chol(sigma));
// }
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::mat sigma) {
   int ncols = sigma.n_cols;
   arma::mat Y = arma::randn(ncols,n);
   return(arma::chol(sigma)*Y);
}

//[[Rcpp::export]]
double get_prior(const arma::vec transformed_theta){
  int num_params = transformed_theta.size();
  // arma::vec mus = {0.2,0.3,10};
  arma::vec transformed_mus = {-0.43275, -5.99021, -1.050226};
  arma::vec transformed_sigmas = {1,1,1};
  double log_prior = 0;
  for(int i=0;i<num_params;++i){
    log_prior += -((transformed_theta[i]-transformed_mus[i])*(transformed_theta[i]-transformed_mus[i]))/(2*transformed_sigmas[i]);
  }
  return(log_prior);
}

//[[Rcpp::export]]
arma::vec transform(const arma::vec &x){
  return(log(exp(x)-1));
}

//[[Rcpp::export]]
arma::vec inv_transform(const arma::vec &x){
  return(log(1+exp(x)));
}