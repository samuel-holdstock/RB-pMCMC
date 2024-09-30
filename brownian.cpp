#include "brownian.h"

//[[Rcpp::export]]
double get_above_maximum_wt(double x,double mu, double wt, double tout){
  if(tout<=1e-15){
    return(0);
  }
  double k = std::max(std::max(wt,0.0),x);
  double probability = std::sqrt(1/(2*M_PI*tout))*std::exp(mu*wt - mu*mu*tout/2 - (2*k-wt)*(2*k-wt)/(2*tout));
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
Rcpp::NumericVector get_expected_jump(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas){
  auto model = model_dict.get_model(model_name);
  Rcpp::NumericMatrix S = model->S;
  Rcpp::NumericVector h = model->get_rates(start,thetas);  
  int num_species = S.nrow();
  int num_reactions = S.ncol();
  int h0 = 0;
  Rcpp::NumericVector expected_jump(num_species);
  for(int i=0;i<num_species;++i){
    for(int j=0;j<num_reactions;++j){
      h0 += h[j];
      expected_jump[i] += 0.5*abs(S(i,j)*h[j])*1;
    }
    expected_jump[i] /= h0;
    h0=0;
  }  
  Rcpp::NumericVector v(num_species);
  for(int i=0;i<num_species;++i){
    // v[i] = 1.28;
    v[i] = 1.28;
  }
  return(v);
  return(expected_jump);
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
  delta = delta+1;
  Rcpp::NumericVector state;
  int state_component;
  Rcpp::NumericVector expected_jump = get_expected_jump(model_name,start,thetas);
  for(int i=0;i<number_species;++i){
    p *= get_probability_hit(start[i],target[i],tout,mu[i],sig[i]);
  }
  for(int i=0;i<total_points;++i){
    state = index_to_state(i,lower,upper);
    for(int j=0;j<number_species;++j){
      a[i] *= get_probability_hit(start[j],state[j],tout-tau,mu[j],sig[j]);
      // std::cout<<"start:"<<start[j]<<", state:"<<state[j]<<std::endl;
      b[i] *= get_probability_hit(state[j],target[j],tau,mu[j],sig[j]);
      alpha[i] = get_above_maximum_sigma_wt(upper[j]+expected_jump[j]-state[j],mu[j],sig[j],target[j]-state[j],tau);
      beta[i] = get_below_minimum_sigma_wt(lower[j]-expected_jump[j]-state[j],mu[j],sig[j],target[j]-state[j],tau);
      delta[i] *= (alpha[i] + beta[i]);
    }
    // std::cout<<"a:"<<a<<std::endl;
    // std::cout<<"b:"<<b<<std::endl;
    // std::cout<<"d:"<<delta<<std::endl;
    // std::cout<<"p:"<<p<<std::endl;
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
  Rcpp::NumericVector expected_jump = get_expected_jump(model_name,start,thetas);
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
    k = std::max(upper[i]+expected_jump[i],target[i]);
    v = std::min(lower[i]-expected_jump[i],target[i]);
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
double get_Q1_brownian_fast(const double &expected_jump, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const double &start, const double &target, 
const double &mu, const double &var){
  double s = tout - tau;
  double pvr1;
  double Q1 = 1;
  pvr1 = R::dnorm(target,start+mu*tout,sqrt(var*tout),0) - 1/sqrt(4*M_PI*var*tau)*R::dnorm(target,start+mu*tout,sqrt(var*(tout-tau/2)),0);
  Q1 *= pvr1;
  return(Q1);
}

//[[Rcpp::export]]
double get_Q2_brownian_fast(const double &expected_jump, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const double &start, const double &target,
const double &upper, const double &mu, const double &var){
  double s = tout - tau;
  double k;
  double pvr2,pvr3,pvr4,pvr5;
  double Q2_a,Q2_b;
  double Q2 = 1;
  k = std::max(upper+expected_jump,target);
  pvr2 = R::dnorm(2*k-target,start+mu*tout,sqrt(var*tout),0)*exp(-2*(k-target)*mu/var)*R::pnorm(k*tout-start*tau,s*(2*k-target),sqrt(tout*var*s*tau),1,0);
  pvr3 = R::dnorm(target,start+mu*tout,sqrt(var*tout),0)*(1-R::pnorm((k*tout-start*tau-target*s)/sqrt(tout*var*s*tau),0,1,1,0));
  pvr4 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(2*k-target,start+mu*tout,sqrt(var*(s+tau/2)),0)*exp(-2*(k-target)*mu/var)*R::pnorm(k*(s+tau/2),(start+mu*s)*tau/2 + (2*k-target-mu*tau)*s,sqrt(var*s*tau/2*(s+tau/2)),1,0);
  pvr5 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(target,start+mu*tout,sqrt(var*(s+tau/2)),0)*(1-R::pnorm(k*(s+tau/2),(start+mu*s)*tau/2+(target-mu*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0));
  Q2_a = pvr2 + pvr3;
  Q2_b = pvr4 + pvr5;
  Q2 *= (Q2_a - Q2_b);
  return(Q2);
}

//[[Rcpp::export]]
double get_Q3_brownian_fast(const double &expected_jump, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const double &start, const double &target, 
const double &lower, const double &mu, const double &var){
  double s = tout - tau;
  double v;
  double pvr6,pvr7,pvr8,pvr9;
  double Q3_a,Q3_b;
  double Q3 = 1;
  v = std::min(lower-expected_jump,target);
  pvr6 = R::dnorm(2*v-target,start+mu*tout,sqrt(var*tout),0)*exp(-2*(v-target)*mu/var)*(1-R::pnorm(v*tout-start*tau,s*(2*v-target),sqrt(tout*var*s*tau),1,0));
  pvr7 = R::dnorm(target,start+mu*tout,sqrt(var*tout),0)*R::pnorm(v*tout-start*tau-target*s,0,sqrt(tout*var*s*tau),1,0);
  pvr8 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(2*v-target,start+mu*tout,sqrt(var*(s+tau/2)),0)*exp(-2*(v-target)*mu/var)*(1-R::pnorm(v*(s+tau/2),(start+mu*s)*tau/2+(2*v-target-mu*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0));
  pvr9 = 1/sqrt(4*M_PI*var*tau)*R::dnorm(target,start+mu*tout,sqrt(var*(s+tau/2)),0)*R::pnorm(v*(s+tau/2),(start+mu*s)*tau/2+(target-mu*tau)*s,sqrt((var*s*tau/2)*(s+tau/2)),1,0);
  Q3_a = pvr6 + pvr7;
  Q3_b = pvr8 + pvr9;
  Q3 *= (Q3_a - Q3_b);
  return(Q3);
}

//[[Rcpp::export]]
Rcpp::List get_box_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target,
const double &pvr_goal){
  Rcpp::NumericVector mu = get_mu(model_name,start,thetas);
  Rcpp::NumericVector var = get_var(model_name,start,thetas);
  Rcpp::NumericVector lower = clone(target);
  Rcpp::NumericVector upper = clone(target);

  int number_species = start.size();
  Rcpp::NumericVector p(number_species);
  Rcpp::NumericVector Q1(number_species);
  Rcpp::NumericVector Q2(number_species);
  Rcpp::NumericVector Q3(number_species);
  Rcpp::NumericVector expected_jump = get_expected_jump(model_name,start,thetas);
  int counter;
  for(int i=0;i<number_species;++i){
    // p[i] = get_probability_hit(start[i],target[i],tout,mu[i],sqrt(var[i]));
    p[i] = 1/sqrt(2*M_PI*var[i]*tout)*exp(-0.5*(target[i]-(start[i]+mu[i]*tout))*(target[i]-(start[i]+mu[i]*tout))/(var[i]*tout)); 
    // Q1[i] = get_Q1_brownian_fast(expected_jump[i],thetas,tout,tau,start[i],target[i],mu[i],var[i]);
    Q1[i] = get_Q1_brownian_fast(0.0,thetas,tout,tau,start[i],target[i],mu[i],var[i]);
    // double goal = Q1[i]*(1 - 1/(2-pvr_goal));
    // double goal = Q1[i]*(1 - pvr_goal)/2;
    double goal = (p[i]*(1-p[i])-Q1[i])*(1 - pvr_goal)/(2*pvr_goal);
    
    counter = 0;
    while(true){
      Q2[i] = get_Q2_brownian_fast(expected_jump[i],thetas,tout,tau,start[i],target[i],upper[i],mu[i],var[i]);
      if(goal>=Q2[i]){
        break;
      }
      upper[i] += 1;
      if((counter+=1)>100){
        std::cout<<"WARNING: BIG BOX, theta:"<<thetas<<", "<<", x:"<<start[i]<<", target:"<<target[i]<<"\n";
        break;
      }
    }
    counter = 0;
    while(true){
        Q3[i] = get_Q3_brownian_fast(expected_jump[i],thetas,tout,tau,start[i],target[i],lower[i],mu[i],var[i]);
      if(goal>=Q3[i]){
        break;
      }
      lower[i] -= 1;
      if((counter+=1)>100){
        std::cout<<"WARNING: BIG BOX, theta:"<<thetas<<", "<<", x:"<<start[i]<<", target:"<<target[i]<<"\n";
        break;
      }
    }

  }
  double Q1_=1;
  double Q2_=1;
  double Q3_=1;
  double p_=1;
  for(int i=0; i<number_species;++i){
    p_ *= p[i];
    Q1_ *= Q1[i];
    Q2_ *= Q2[i];
    Q3_ *= Q3[i];
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("lower")=lower,Rcpp::Named("upper")=upper,Rcpp::Named("estimated_PVR_inf_box")=(Q1)/(p_*(1-p_)));
  return(results);
}

//[[Rcpp::export]]
Rcpp::List get_box_brownian_fast_list(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &tout_list, const Rcpp::NumericVector &tau_list, const Rcpp::NumericVector &start, const Rcpp::NumericMatrix &target_list,
const double &pvr_goal){
  int num_obs = tout_list.size();
  int num_species = start.size();
  double tout;
  double tau;
  double prev_tout=0;
  Rcpp::NumericVector x = start;
  Rcpp::NumericVector obs(num_species);
  Rcpp::NumericMatrix lower_list(num_obs,num_species);
  Rcpp::NumericMatrix upper_list(num_obs,num_species);
  Rcpp::NumericVector lower(num_species);
  Rcpp::NumericVector upper(num_species);
  Rcpp::DataFrame box;
  for(int i=0;i<num_obs;++i){
    tout = tout_list[i]-prev_tout;
    tau = tau_list[i];    
    obs = target_list(i,Rcpp::_);
    box = get_box_brownian_fast(model_name,thetas,tout,tau,x,obs,pvr_goal);
    lower = box["lower"];
    upper = box["upper"];
    lower_list(i,Rcpp::_) = lower;
    upper_list(i,Rcpp::_) = upper;
    x = clone(obs);
    prev_tout = tout_list[i];
  }
  Rcpp::List list = Rcpp::List::create(Rcpp::Named("lower")=lower_list, Rcpp::Named("upper")=upper_list );
  // Rcpp::List list = Rcpp::List::create(Rcpp::Named("lower")=0, Rcpp::Named("upper")=0 );
  return(list);
}

//[[Rcpp::export]]
double get_var_big_box_tau(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &target, const double &tout, const double tau){
  double mu = get_mu(model_name,start,thetas)[0];
  double sig2 = get_var(model_name,start,thetas)[0];
  double RB_variance = 1/sqrt(4*M_PI*sig2*tau)*1/sqrt(2*M_PI*sig2*(tout-tau*0.5))*exp(-0.5*(target[0]-(start[0]+mu*tout))*(target[0]-(start[0]+mu*tout))/(sig2*(tout-tau*0.5)))-
    1/(2*M_PI*sig2*tout)*exp(-0.5*(target[0]-start[0]-mu*tout)*(target[0]-start[0]-mu*tout)/(sig2*tout/2));
  double p = 1/sqrt(2*M_PI*sig2*tout)*exp(-0.5*(target[0]-start[0]-mu*tout)*(target[0]-start[0]-mu*tout)/(sig2*tout));
  if(RB_variance>p*(1-p)){
    RB_variance = p*(1-p);
  }
  return(RB_variance);
}

//[[Rcpp::export]]
double get_PVR_big_box_tau(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &target, const double &tout, const double tau){
  double mu = get_mu(model_name,start,thetas)[0];
  double sig2 = get_var(model_name,start,thetas)[0];
  double RB_variance = 1/sqrt(4*M_PI*sig2*tau)*1/sqrt(2*M_PI*sig2*(tout-tau/2))*exp(-0.5*(target[0]-(start[0]+mu*tout))*(target[0]-(start[0]+mu*tout))/(sig2*(tout-tau/2)))-
    1/(2*M_PI*sig2*tout)*exp(-0.5*(target[0]-start[0]-mu*tout)*(target[0]-start[0]-mu*tout)/(sig2*tout/2));
  double p = 1/sqrt(2*M_PI*sig2*tout)*exp(-0.5*(target[0]-start[0]-mu*tout)*(target[0]-start[0]-mu*tout)/(sig2*tout));
  double PVR = (p*(1-p)-RB_variance)/(p*(1-p));
  if(PVR<0){
    PVR = 0;
  }
  return(PVR);
}

int sign(const double &x){
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

//[[Rcpp::export]]
double get_tau(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericVector &target, const double &tout, const double &goal){
  double a = 0;
  double b = tout;
  double tol=0.00001;
  double tau;
  int counter = 0;
  while (abs(b - a) > tol){
    tau = (a + b) / 2;
    if(get_PVR_big_box_tau(model_name,start,thetas,target,tout,tau) == goal){
      return(tau);
    }
    if(sign(get_PVR_big_box_tau(model_name,start,thetas,target,tout,a)-goal) == sign(get_PVR_big_box_tau(model_name,start,thetas,target,tout,tau)-goal)){
      a = tau;
    }
    else{
      b = tau;
    }
    counter += 1;
    if(counter>1000){
      std::cout<<"WARNING: TAU NOT FOUND. A:"<<a<<",B:"<<b<<std::endl;
      return(-1);
    }
  }
  return(b);
}

//[[Rcpp::export]]
Rcpp::NumericVector get_tau_list(const std::string &model_name, const Rcpp::NumericVector &start, const Rcpp::NumericVector &thetas, 
const Rcpp::NumericMatrix &obs_list, const Rcpp::NumericVector &tout_list, const double &goal){
  int num_obs = obs_list.size();
  Rcpp::NumericVector target(num_obs);
  Rcpp::NumericVector x = start;
  Rcpp::NumericVector tau_list(num_obs);
  double tout;
  double prev_tout=0;
  for(int i=0;i<num_obs;++i){
    target = obs_list(i,Rcpp::_);
    tout = tout_list[i];
    tau_list[i] = get_tau(model_name,x,thetas,target,tout-prev_tout,goal);
    x = clone(target);
    prev_tout=tout;
  }
  return(tau_list);
}


// Rcpp::DataFrame get_box_brownian_fast(const std::string &model_name, const Rcpp::NumericVector &thetas, 
// const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target,
// const double &pvr_goal){
//   Rcpp::NumericVector mu = get_mu(model_name,start,thetas);
//   Rcpp::NumericVector var = get_var(model_name,start,thetas);
//   Rcpp::NumericVector lower = clone(target);
//   Rcpp::NumericVector upper = clone(target);
//   double p = 1;
//   int number_species = start.size();
//   for(int i=0;i<number_species;++i){
//     p *= get_probability_hit(start[i],target[i],tout,mu[i],sqrt(var[i]));
//   }
//   double Q1 = get_Q1_brownian_fast(model_name,thetas,tout,tau,start,target,mu,var);
//   double Q2 = 0;
//   double Q3 = 0;
//   //std::cout<<"Upper:";
//   //double goal = (Q1 - pvr_goal*p*(1-p))/2;
//   double goal = Q1*(1 - pvr_goal)/2;
//   // Q2/Q1 <= (1-pvr_goal)/2
  
//   int counter = 0;
//   while(true){
//     //std::cout<<upper[0]<<", ";
//     Q2 = get_Q2_brownian_fast(model_name,thetas,tout,tau,start,target,upper,mu,var);
//     if(goal>=Q2){
//       break;
//     }
//     upper[0] += 1;
//     if((counter+=1)>100){
//       std::cout<<"WARNING: BIG BOX\n";
//       break;
//     }
//   }
//   counter = 0;
//   //std::cout<<std::endl;
//   //std::cout<<"Lower:";
//   while(true){
//     //std::cout<<lower[0]<<", ";
//     Q3 = get_Q3_brownian_fast(model_name,thetas,tout,tau,start,target,lower,mu,var);
//     if(goal>=Q3){
//       break;
//     }
//     lower[0] -= 1;
//     if((counter+=1)>100){
//       std::cout<<"WARNING: BIG BOX\n";
//       break;
//     }
//   }
//   //std::cout<<std::endl;
//   //std::cout<<"Q1:"<<Q1<<", Q2:"<<Q2<<", Q3:"<<Q3<<std::endl;
//   //Rcpp::List results = Rcpp::List::create(Rcpp::Named("Q1")=Q1,Rcpp::Named("Q2")=Q2,Rcpp::Named("Q3")=Q3);
//   Rcpp::DataFrame results = Rcpp::DataFrame::create(Rcpp::Named("lower")=lower,Rcpp::Named("upper")=upper,Rcpp::Named("scale")=p*(1-p)/(Q1*(1-Q1)));
//   return(results);
// }
