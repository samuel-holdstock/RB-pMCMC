#include "brownian.h"

//[[Rcpp::export]]
double get_above_maximum_wt(double x,double mu, double wt, double tout){
  if(tout<=1e-15){
    return(0);
  }
  x = std::max(std::max(wt,0.0),x);
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
      // std::cout<<"start:"<<start[j]<<", state:"<<state[j]<<std::endl;
      b[i] *= get_probability_hit(state[j],target[j],tau,mu[j],sig[j]);
      alpha[i] = get_above_maximum_sigma_wt(upper[j]+1-state[j],mu[j],sig[j],target[j]-state[j],tau);
      beta[i] = get_below_minimum_sigma_wt(lower[j]-1-state[j],mu[j],sig[j],target[j]-state[j],tau);
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
    // percentage_variance_reduction -= (a[i]*b[i]*(1-b[i]) - (a[i]*(alpha[i]*(1-alpha[i])+beta[i]*(1-beta[i]))));
    
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
Rcpp::DataFrame get_variance_brownian2(const std::string &model_name, const Rcpp::NumericVector &thetas, 
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
      // std::cout<<"start:"<<start[j]<<", state:"<<state[j]<<std::endl;
      b[i] *= get_probability_hit(state[j],target[j],tau,mu[j],sig[j]);
      alpha[i] = get_above_maximum_sigma_wt(upper[j]+1-state[j],mu[j],sig[j],target[j]-state[j],tau);
      beta[i] = get_below_minimum_sigma_wt(lower[j]-1-state[j],mu[j],sig[j],target[j]-state[j],tau);
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
    // percentage_variance_reduction -= (a[i]*b[i]*(1-b[i]) - (a[i]*(alpha[i]*(1-alpha[i])+beta[i]*(1-beta[i]))));
    
  }
  //double pvr = 1 - (Q1 - Q2 - Q3)/(p*(1-p));
  //std::cout<<"PVR:"<<pvr<<", Q1: "<<Q1<<", Q2:"<<Q2<<", Q3:"<<Q3<<", p:"<<p<<std::endl;
  percentage_variance_reduction = percentage_variance_reduction/(p*(1-p));
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
    Rcpp::Named("a")=a,
    Rcpp::Named("b")=b,
    Rcpp::Named("alpha")=alpha,
    Rcpp::Named("beta")=beta,
    Rcpp::Named("delta")=delta,
    Rcpp::Named("p")=p,
    Rcpp::Named("pvr")=percentage_variance_reduction
    );
  return(df);
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
    return(1);
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
  double k,v;
  double var;
  double pvr1,pvr2,pvr3,pvr4,pvr5,pvr6,pvr7,pvr8,pvr9;
  double Q1 = 1;
  double Q2_a,Q2_b;
  double Q2 = 0;
  double Q3_a,Q3_b;
  double Q3 = 0;
  for(int i=0;i<number_species;++i){
    k = std::max(upper[i]+1,target[i]);
    v = std::min(lower[i]-1,target[i]);
    var = sig[i]*sig[i];
    pvr1 = R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var*tout),0) - 1/sqrt(4*M_PI*var*tau)*R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var*(s+tau/2)),0);
    pvr2 = R::dnorm(2*k-target[i],start[i]+mu[i]*tout,sqrt(var*tout),0)*exp(-2*(k-target[i])*mu[i]/var)*R::pnorm(k*tout-start[i]*tau,s*(2*k-target[i]),sqrt(tout*var*s*tau),1,0);
    pvr3 = R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var*tout),0)*(1-R::pnorm((k*tout-start[i]*tau-target[i]*s)/sqrt(tout*var*s*tau),0,1,1,0));
    pvr4 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(2*k-target[i],start[i]+mu[i]*tout,sqrt(var*(s+tau/2)),0)*exp(-2*(k-target[i])*mu[i]/var)*R::pnorm(k*(s+tau/2),(start[i]+mu[i]*s)*tau/2 + (2*k-target[i]-mu[i]*tau)*s,sqrt(var*s*tau/2*(s+tau/2)),1,0);
    pvr5 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var*(s+tau/2)),0)*(1-R::pnorm(k*(s+tau/2),(start[i]+mu[i]*s)*tau/2+(target[i]-mu[i]*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0));
    pvr6 = R::dnorm(2*v-target[i],start[i]+mu[i]*tout,sqrt(var*tout),0)*exp(-2*(v-target[i])*mu[i]/var)*(1-R::pnorm(v*tout-start[i]*tau,s*(2*v-target[i]),sqrt(tout*var*s*tau),1,0));
    pvr7 = R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var*tout),0)*R::pnorm(v*tout-start[i]*tau-target[i]*s,0,sqrt(tout*var*s*tau),1,0);
    pvr8 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(2*v-target[i],start[i]+mu[i]*tout,sqrt(var*(s+tau/2)),0)*exp(-2*(v-target[i])*mu[i]/var)*(1-R::pnorm(v*(s+tau/2),(start[i]+mu[i]*s)*tau/2+(2*v-target[i]-mu[i]*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0));
    pvr9 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var*(s+tau/2)),0)*R::pnorm(v*(s+tau/2),(start[i]+mu[i]*s)*tau/2+(target[i]-mu[i]*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0);
    Q1 *= pvr1;
    Q2_a = pvr2 + pvr3;
    Q2_b = pvr4 + pvr5;
    Q3_a = pvr6 + pvr7;
    Q3_b = pvr8 + pvr9;
    Q2 += Q2_a - Q2_b;
    Q3 += Q3_a - Q3_b;
  }
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
//[[Rcpp::export]]
double get_Q1_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &mu, const Rcpp::NumericVector &var){
  double s = tout - tau;
  double pvr1;
  double Q1 = 1;
  int number_species = start.size();
  for(int i=0;i<number_species;++i){
    pvr1 = R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var[i]*tout),0) - 1/sqrt(4*M_PI*var[i]*tau)*R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var[i]*(s+tau/2)),0);
    Q1 *= pvr1;
  }
  return(Q1);
}

//[[Rcpp::export]]
double get_Q2_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target,
const Rcpp::NumericVector &upper, const Rcpp::NumericVector &mu, const Rcpp::NumericVector &var){
  double s = tout - tau;
  double k;
  double pvr2,pvr3,pvr4,pvr5;
  double Q2_a,Q2_b;
  double Q2 = 0;
  int number_species = start.size();
  for(int i=0;i<number_species;++i){
    k = std::max(upper[i],target[i]);
    pvr2 = R::dnorm(2*k-target[i],start[i]+mu[i]*tout,sqrt(var[i]*tout),0)*exp(-2*(k-target[i])*mu[i]/var[i])*R::pnorm(k*tout-start[i]*tau,s*(2*k-target[i]),sqrt(tout*var[i]*s*tau),1,0);
    pvr3 = R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var[i]*tout),0)*(1-R::pnorm((k*tout-start[i]*tau-target[i]*s)/sqrt(tout*var[i]*s*tau),0,1,1,0));
    pvr4 = 1/sqrt(4*M_PI*var[i]*tau)*R::dnorm(2*k-target[i],start[i]+mu[i]*tout,sqrt(var[i]*(s+tau/2)),0)*exp(-2*(k-target[i])*mu[i]/var[i])*R::pnorm(k*(s+tau/2),(start[i]+mu[i]*s)*tau/2 + (2*k-target[i]-mu[i]*tau)*s,sqrt(var[i]*s*tau/2*(s+tau/2)),1,0);
    pvr5 = 1/sqrt(4*M_PI*var[i]*tau)*R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var[i]*(s+tau/2)),0)*(1-R::pnorm(k*(s+tau/2),(start[i]+mu[i]*s)*tau/2+(target[i]-mu[i]*tau)*s,sqrt((var[i]*s*tau/2)*(s+tau/2)),1,0));
    Q2_a = pvr2 + pvr3;
    Q2_b = pvr4 + pvr5;
    Q2 += Q2_a - Q2_b;
  }
  return(Q2);
}

//[[Rcpp::export]]
double get_Q3_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &lower, const Rcpp::NumericVector &mu, const Rcpp::NumericVector &var){
  double s = tout - tau;
  double v;
  double pvr6,pvr7,pvr8,pvr9;
  double Q3_a,Q3_b;
  double Q3 = 0;
  int number_species = start.size();
  for(int i=0;i<number_species;++i){
    v = std::min(lower[i],target[i]);
    pvr6 = R::dnorm(2*v-target[i],start[i]+mu[i]*tout,sqrt(var[i]*tout),0)*exp(-2*(v-target[i])*mu[i]/var[i])*(1-R::pnorm(v*tout-start[i]*tau,s*(2*v-target[i]),sqrt(tout*var[i]*s*tau),1,0));
    pvr7 = R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var[i]*tout),0)*R::pnorm(v*tout-start[i]*tau-target[i]*s,0,sqrt(tout*var[i]*s*tau),1,0);
    pvr8 = 1/sqrt(4*M_PI*var[i]*tau)*R::dnorm(2*v-target[i],start[i]+mu[i]*tout,sqrt(var[i]*(s+tau/2)),0)*exp(-2*(v-target[i])*mu[i]/var[i])*(1-R::pnorm(v*(s+tau/2),(start[i]+mu[i]*s)*tau/2+(2*v-target[i]-mu[i]*tau)*s,sqrt((var[i]*s*tau/2)*(s+tau/2)),1,0));
    pvr9 = 1/sqrt(4*M_PI*var[i]*tau)*R::dnorm(target[i],start[i]+mu[i]*tout,sqrt(var[i]*(s+tau/2)),0)*R::pnorm(v*(s+tau/2),(start[i]+mu[i]*s)*tau/2+(target[i]-mu[i]*tau)*s,sqrt((var[i]*s*tau/2)*(s+tau/2)),1,0);
    Q3_a = pvr6 + pvr7;
    Q3_b = pvr8 + pvr9;
    Q3 += Q3_a - Q3_b;
  }
  return(Q3);
}

//[[Rcpp::export]]
Rcpp::DataFrame get_box_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target,
const double &pvr_goal){
  Rcpp::NumericVector mu = get_mu(model_name,start,thetas);
  Rcpp::NumericVector var = get_var(model_name,start,thetas);
  Rcpp::NumericVector lower = clone(target);
  Rcpp::NumericVector upper = clone(target);
  double p = 1;
  int number_species = start.size();
  for(int i=0;i<number_species;++i){
    p *= get_probability_hit(start[i],target[i],tout,mu[i],sqrt(var[i]));
  }
  double Q1 = get_Q1_brownian_fast(model_name,thetas,tout,tau,start,target,mu,var);
  double Q2 = 0;
  double Q3 = 0;
  //std::cout<<"Upper:";
  //double goal = (Q1 - pvr_goal*p*(1-p))/2;
  double goal = Q1*(1 - pvr_goal)/2;
  // Q2/Q1 <= (1-pvr_goal)/2
  
  int counter = 0;
  while(true){
    //std::cout<<upper[0]<<", ";
    Q2 = get_Q2_brownian_fast(model_name,thetas,tout,tau,start,target,upper,mu,var);
    if(goal>=Q2){
      break;
    }
    upper[0] += 1;
    if((counter+=1)>100){
      std::cout<<"ERROR";
      break;
    }
  }
  counter = 0;
  //std::cout<<std::endl;
  //std::cout<<"Lower:";
  while(true){
    //std::cout<<lower[0]<<", ";
    Q3 = get_Q3_brownian_fast(model_name,thetas,tout,tau,start,target,lower,mu,var);
    if(goal>=Q3){
      break;
    }
    lower[0] -= 1;
    if((counter+=1)>100){
      std::cout<<"ERROR";
      break;
    }
  }
  //std::cout<<std::endl;
  //std::cout<<"Q1:"<<Q1<<", Q2:"<<Q2<<", Q3:"<<Q3<<std::endl;
  //Rcpp::List results = Rcpp::List::create(Rcpp::Named("Q1")=Q1,Rcpp::Named("Q2")=Q2,Rcpp::Named("Q3")=Q3);
  Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("lower")=lower,Rcpp::Named("upper")=upper);
  return(results);
}
