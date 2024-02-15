#include "exact.h"

//[[Rcpp::export]]
double get_variance_exact(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, 
const Rcpp::NumericVector &lower_limit, const Rcpp::NumericVector &upper_limit){
  Rcpp::NumericMatrix Q = get_coffin_matrix(model_name,lower,upper,thetas);
  Rcpp::NumericMatrix Qt = Rcpp::transpose(Q);  
  Rcpp::NumericMatrix large_Q = get_coffin_matrix(model_name,lower_limit,upper_limit,thetas);
  Rcpp::NumericMatrix large_Qt = Rcpp::transpose(large_Q);
  arma::mat v(Qt.nrow(),1);
  arma::mat w(1,Q.nrow());
  arma::mat x(large_Qt.nrow(),1);
  arma::mat y(1,large_Q.nrow());
  int obs_index = state_to_index(target,lower,upper);
  int start_index = state_to_index(start,lower,upper);
  v[obs_index] = 1;
  arma::mat delta = vT_exp_Q(v,Qt*tau,1e-20,false,true,false);  // actually Q

  start_index = state_to_index(start,lower_limit,upper_limit);
  obs_index = state_to_index(target,lower_limit,upper_limit);
  y[start_index] = 1;
  x[obs_index] = 1;
  arma::mat a = v_exp_Q(y,large_Q*(tout-tau),1e-20,false,true,false); 
  arma::mat b = vT_exp_Q(x,large_Qt*tau,1e-20,false,true,false); 
  double p = (v_exp_Q(y,large_Q*tout,1e-20,false,true,false))[obs_index];

  double percentage_variance_reduction = p*(1-p);
  Rcpp::NumericVector a_(Q.nrow()-1);
  Rcpp::NumericVector b_(Q.nrow()-1);
  int counter = 0;
  for(int i=0;i<large_Q.nrow()-1;++i){
    if(!(inBox(index_to_state(i,lower_limit,upper_limit),lower,upper))){
      continue;
    }
    a_[counter] = a[i];
    b_[counter] = b[i];
    //std::cout<<"counter:"<<index_to_state(counter,lower,upper)<<" index:"<<index_to_state(i,lower_limit,upper_limit)<<std::endl;
    //std::cout<<"a:"<<a[i]<<", b:"<<b[i]<<", delta"<<b[i]-delta[counter]<<std::endl;
    percentage_variance_reduction -= (a[i]*b[i]*(1-b[i])-a[i]*(b[i]-delta[counter])*(1-(b[i]-delta[counter]))); 
    counter+=1;
  }
  percentage_variance_reduction /= p*(1-p);
  return(percentage_variance_reduction);
}
//[[Rcpp::export]]
Rcpp::DataFrame get_variance_exact2(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const double &tau, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, 
const Rcpp::NumericVector &lower_limit, const Rcpp::NumericVector &upper_limit){
  Rcpp::NumericMatrix Q = get_coffin_matrix(model_name,lower,upper,thetas);
  Rcpp::NumericMatrix Qt = Rcpp::transpose(Q);  
  Rcpp::NumericMatrix large_Q = get_coffin_matrix(model_name,lower_limit,upper_limit,thetas);
  Rcpp::NumericMatrix large_Qt = Rcpp::transpose(large_Q);
  arma::mat v(Qt.nrow(),1);
  arma::mat w(1,Q.nrow());
  arma::mat x(large_Qt.nrow(),1);
  arma::mat y(1,large_Q.nrow());
  int obs_index = state_to_index(target,lower,upper);
  int start_index = state_to_index(start,lower,upper);
  v[obs_index] = 1;
  arma::mat delta = vT_exp_Q(v,Qt*tau,1e-20,false,true,false);  // actually Q

  start_index = state_to_index(start,lower_limit,upper_limit);
  obs_index = state_to_index(target,lower_limit,upper_limit);
  y[start_index] = 1;
  x[obs_index] = 1;
  arma::mat a = v_exp_Q(y,large_Q*(tout-tau),1e-20,false,true,false); 
  arma::mat b = vT_exp_Q(x,large_Qt*tau,1e-20,false,true,false); 
  double p = (v_exp_Q(y,large_Q*tout,1e-20,false,true,false))[obs_index];

  double percentage_variance_reduction = p*(1-p);
  Rcpp::NumericVector a2(Q.nrow()-1);
  Rcpp::NumericVector b2(Q.nrow()-1);
  Rcpp::NumericVector delta2(Q.nrow()-1);
  int counter = 0;
  for(int i=0;i<large_Q.nrow()-1;++i){
    if(!(inBox(index_to_state(i,lower_limit,upper_limit),lower,upper))){
      continue;
    }
    a2[counter] = a[i];
    b2[counter] = b[i];
    //std::cout<<"counter:"<<index_to_state(counter,lower,upper)<<" index:"<<index_to_state(i,lower_limit,upper_limit)<<std::endl;
    std::cout<<index_to_state(counter,lower,upper)<<", a:"<<a[i]<<", b:"<<b[i]<<", delta"<<b[i]-delta[counter]<<std::endl;
    percentage_variance_reduction -= (a[i]*b[i]*(1-b[i])-a[i]*(b[i]-delta[counter])*(1-(b[i]-delta[counter]))); 
    delta2[counter] = b[i]-delta[counter];
    if(delta2[counter]<0){
      delta2[counter] = 0;
    }
    //std::cout<<index_to_state(counter,lower,upper)<<", "<<delta2[counter]<<std::endl;
    counter+=1;
  }
  percentage_variance_reduction /= p*(1-p);
  Rcpp::DataFrame df = Rcpp::DataFrame::create(
    Rcpp::Named("a")=a2,
    Rcpp::Named("b")=b2,
    Rcpp::Named("delta")=delta2,
    Rcpp::Named("p")=p,
    Rcpp::Named("pvr")=percentage_variance_reduction
  );
  return(df);
  return(percentage_variance_reduction);
}


//[[Rcpp::export]]
Rcpp::NumericMatrix get_contour_exact(const std::string &model_name, const Rcpp::NumericVector &thetas, 
const double &tout, const Rcpp::NumericVector &taus, const Rcpp::NumericVector &start, const Rcpp::NumericVector &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list, 
const Rcpp::NumericVector &lower_limit, const Rcpp::NumericVector &upper_limit){
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
      z(j,i) = 1-get_variance_exact(model_name, thetas, tout, tau, start, target, lower, upper,lower_limit,upper_limit);
    }
  }
  return(z);
}

