#include "brownian.h"

//[[Rcpp::export]]
double get_above_maximum_wt(double x,double mu, double wt, double tout){
  if(tout<=1e-15){
    return(0);
  }
  double probability = std::sqrt(1/(2*M_PI*tout))*std::exp(mu*wt - std::pow(mu,2)*tout/2 - (std::pow(2*x-wt,2))/(2*tout));
  //std::cout<<"x:"<<x<<", mu:"<<mu<<", wt:"<<wt<<", tout:"<<tout<<std::endl;
  return(probability);
}

//[[Rcpp::export]]
double get_above_maximum_sigma_wt(double x, double mu, double sig, double wt, double tout){
  double probability = get_above_maximum_wt(x/sig,mu/sig,wt/sig,tout)/sig;
  return(probability);
}

//[[Rcpp::export]]
double get_below_minimum_sigma_wt(double x, double mu, double sig, double wt, double tout){ 
  double probability = get_above_maximum_wt(-x/sig,-mu/sig,-wt/sig,tout)/sig;
  return(probability);
}

//[[Rcpp::export]]
double get_probability_hit(double start, double target, double tout, double mu, double sig){
  double probability = R::pnorm(target + 0.5, start + mu*tout, sqrt(tout)*sig,true,false) -
    R::pnorm(target - 0.5, start + mu*tout, sqrt(tout)*sig,true,false);
  return(probability);
}

//[[Rcpp::export]]
Rcpp::NumericVector get_mu(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas){
  auto model = model_dict.get_model(model_name);
  Rcpp::NumericMatrix S = model->S;
  Rcpp::NumericVector h = model->get_rates(start,thetas);  
  int number_reactions = S.ncol();
  int number_species = S.nrow();
  Rcpp::NumericVector mu(number_species);
  for(int i=0;i<number_species;++i){
    for(int j=0;j<number_reactions;++j){
      mu(i) += S(i,j)*h[j];
    }
  }  
  return(mu);
}

//[[Rcpp::export]]
Rcpp::NumericVector get_var(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas){
  auto model = model_dict.get_model(model_name);
  Rcpp::NumericMatrix S = model->S;
  Rcpp::NumericVector h = model->get_rates(start,thetas);  
  int number_reactions = S.ncol();
  int number_species = S.nrow();
  Rcpp::NumericVector sig(number_species);
  for(int i=0;i<number_species;++i){
    for(int j=0;j<number_reactions;++j){
      sig(i) += S(i,j)*h[j]*S(i,j);
      //sig(i,i) += S(i,j)*h[j]*S(i,j);
    }
  }  
  return(sig);
}

//[[Rcpp::export]]
double get_variance_brownian(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper){
  auto model = model_dict.get_model(model_name);
  Rcpp::NumericMatrix S = model->S;
  Rcpp::NumericVector h = model->get_rates(start,thetas);  
  int number_reactions = S.ncol();
  int number_species = S.nrow();
  Rcpp::NumericVector mu(number_species);
  Rcpp::NumericVector sig(number_species);
  for(int i=0;i<number_species;++i){
    for(int j=0;j<number_reactions;++j){
      mu(i) += S(i,j)*h[j];
      sig(i) += S(i,j)*h[j]*S(i,j);
    }
    sig(i) = sqrt(sig(i));
  }  
  int total_points = 1;
  for(int i=0;i<number_species;++i){
      total_points *= (upper[i]-lower[i]+1);
  }
  Rcpp::NumericVector a(total_points);
  Rcpp::NumericVector b(total_points);
  Rcpp::NumericVector delta(total_points);
  Rcpp::NumericVector alpha(total_points);
  Rcpp::NumericVector beta(total_points);
  double p = 1;
  a = a+1;
  b = b+1;
  Rcpp::NumericVector state;
  int state_component;
  for(int i=0;i<number_species;++i){
    p *= get_probability_hit(start[i],target[i],tout,mu[i],sig[i]);
  }
  for(int i=0;i<total_points;++i){
    state = index_to_state(i,lower,upper);
    for(int j=0;j<number_species;++j){
      a[i] *= get_probability_hit(start[j],state[j],tout-tau,mu[j],sig[j]);
      b[i] *= get_probability_hit(state[j],target[j],tau,mu[j],sig[j]);
      alpha[i] = get_above_maximum_sigma_wt(upper[j]-state[j]+1,mu[j],sig[j],target[j]-state[j],tau);
      beta[i] = get_below_minimum_sigma_wt(-(state[j]-lower[j]+1),mu[j],sig[j],target[j]-state[j],tau);
      delta[i] += (alpha[i] + beta[i]);
    }
  }

  double percentage_variance_reduction = p*(1-p);
  // double Q1 = 0;
  // double Q2 = 0;
  // double Q3 = 0;
  for(int i=0;i<total_points;++i){
    // Q1 += a[i]*b[i]*(1-b[i]);
    // Q2 += a[i]*alpha[i]*(1-alpha[i]);
    // Q3 += a[i]*beta[i]*(1-beta[i]);
    // std::cout<<"a:"<<a[i]<<", b:"<<b[i]<<", delta:"<<delta[i]<<std::endl;
    // std::cout<<"alpha:"<<alpha[i]<<", beta: "<<beta[i]<<std::endl;
    percentage_variance_reduction -= (a[i]*b[i]*(1-b[i]) - a[i]*delta[i]*(1-delta[i]));
  }
  //double pvr = 1 - (Q1 - Q2 - Q3)/(p*(1-p));
  //std::cout<<"PVR:"<<pvr<<", Q1: "<<Q1<<", Q2:"<<Q2<<", Q3:"<<Q3<<", p:"<<p<<std::endl;
  percentage_variance_reduction = percentage_variance_reduction/(p*(1-p));
  if(percentage_variance_reduction>1){
    return(1);
  }
  if(percentage_variance_reduction<0){
    return(0);
  }
  return(percentage_variance_reduction);
}

//[[Rcpp::export]]
Rcpp::NumericMatrix get_contour_brownian(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const Rcpp::NumericVector &taus, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list){
  int number_widths = lower_list.size();
  int number_tau = taus.size();
  Rcpp::NumericVector lower;
  Rcpp::NumericVector upper;
  double tau;
  Rcpp::NumericMatrix z(number_tau,number_widths);  
  for(int i=0;i<number_widths;++i){
    for(int j=0;j<number_tau;++j){
      lower = lower_list[i];
      upper = upper_list[i];
      tau = taus[j];
      z(j,i) = 1-get_variance_brownian(model_name, thetas, tout, tau, start, target, lower, upper);
    }
  }
  return(z);
}

//[[Rcpp::export]]
double get_variance_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper){
  if(tau<1e-15){
    return(0);
  }
  auto model = model_dict.get_model(model_name);
  Rcpp::NumericMatrix S = model->S;
  Rcpp::NumericVector h = model->get_rates(start,thetas);  
  int number_reactions = S.ncol();
  int number_species = S.nrow();
  Rcpp::NumericVector mu(number_species);
  Rcpp::NumericVector sig(number_species);
  for(int i=0;i<number_species;++i){
    for(int j=0;j<number_reactions;++j){
      mu(i) += S(i,j)*h[j];
      sig(i) += S(i,j)*h[j]*S(i,j);
    }
    sig(i) = sqrt(sig(i));
  }  
  int total_points = 1;
  double p = 1;
  for(int i=0;i<number_species;++i){
    total_points *= (upper[i]-lower[i]+1);
    p *= get_probability_hit(start[i],target[i],tout,mu[i],sig[i]);
  }
  double percentage_variance_reduction;
  double s = tout - tau;
  double k = std::max(upper[0],target[0]);
  double v = std::min(lower[0],target[0]);
  double var = sig[0]*sig[0];
  double pvr1 = R::dnorm(target[0],start[0]+mu[0]*tout,sqrt(var*tout),0) - 1/sqrt(4*M_PI*var*tau)*R::dnorm(target[0],start[0]+mu[0]*tout,sqrt(var*(s+tau/2)),0);
  double pvr2 = R::dnorm(2*k-target[0],start[0]+mu[0]*tout,sqrt(var*tout),0)*exp(-2*(k-target[0])*mu[0]/var)*R::pnorm(k*tout-start[0]*tau,s*(2*k-target[0]),sqrt(tout*var*s*tau),1,0);
  double pvr3 = R::dnorm(target[0],start[0]+mu[0]*tout,sqrt(var*tout),0)*(1-R::pnorm((k*tout-start[0]*tau-target[0]*s)/sqrt(tout*var*s*tau),0,1,1,0));
  double pvr4 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(2*k-target[0],start[0]+mu[0]*tout,sqrt(var*(s+tau/2)),0)*exp(-2*(k-target[0])*mu[0]/var)*R::pnorm(k*(s+tau/2),(start[0]+mu[0]*s)*tau/2 + (2*k-target[0]-mu[0]*tau)*s,sqrt(var*s*tau/2*(s+tau/2)),1,0);
  double pvr5 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(target[0],start[0]+mu[0]*tout,sqrt(var*(s+tau/2)),0)*(1-R::pnorm(k*(s+tau/2),(start[0]+mu[0]*s)*tau/2+(target[0]-mu[0]*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0));
  double pvr6 = R::dnorm(2*v-target[0],start[0]+mu[0]*tout,sqrt(var*tout),0)*exp(-2*(v-target[0])*mu[0]/var)*(1-R::pnorm(v*tout-start[0]*tau,s*(2*v-target[0]),sqrt(tout*var*s*tau),1,0));
  double pvr7 = R::dnorm(target[0],start[0]+mu[0]*tout,sqrt(var*tout),0)*R::pnorm(v*tout-start[0]*tau-target[0]*s,0,sqrt(tout*var*s*tau),1,0);
  double pvr8 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(2*v-target[0],start[0]+mu[0]*tout,sqrt(var*(s+tau/2)),0)*exp(-2*(v-target[0])*mu[0]/var)*(1-R::pnorm(v*(s+tau/2),(start[0]+mu[0]*s)*tau/2+(2*v-target[0]-mu[0]*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0));
  double pvr9 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(target[0],start[0]+mu[0]*tout,sqrt(var*(s+tau/2)),0)*R::pnorm(v*(s+tau/2),(start[0]+mu[0]*s)*tau/2+(target[0]-mu[0]*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0);
  double Q1 = pvr1;
  double Q2_a = pvr2 + pvr3;
  double Q2_b = pvr4 + pvr5;
  double Q3_a = pvr6 + pvr7;
  double Q3_b = pvr8 + pvr9;
  double Q2 = Q2_a - Q2_b;
  double Q3 = Q3_a - Q3_b;
  //std::cout<<"a:"<<a_<<", b:"<<b_<<std::endl;
  percentage_variance_reduction = 1 - (Q1 - Q2 - Q3)/(p*(1-p));
  //std::cout<<"pvr:"<<percentage_variance_reduction<<", Q1:"<<Q1<<", Q2:"<<Q2<<", Q3:"<<Q3<<", p:"<<p<<std::endl;
  if(percentage_variance_reduction>1){
    percentage_variance_reduction = 1;
  }
  if(percentage_variance_reduction<0){
    percentage_variance_reduction = 0;
  }
  return(percentage_variance_reduction);
}

//[[Rcpp::export]]
Rcpp::NumericMatrix get_contour_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const Rcpp::NumericVector &taus, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list){
  int number_widths = lower_list.size();
  int number_tau = taus.size();
  Rcpp::NumericVector lower;
  Rcpp::NumericVector upper;
  double tau;
  Rcpp::NumericMatrix z(number_tau,number_widths);  
  for(int i=0;i<number_widths;++i){
    for(int j=0;j<number_tau;++j){
      lower = lower_list[i];
      upper = upper_list[i];
      tau = taus[j];
      z(j,i) = 1-get_variance_brownian_fast(model_name, thetas, tout, tau, start, target, lower, upper);
    }
  }
  return(z);
}

