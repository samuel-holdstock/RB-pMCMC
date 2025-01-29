#include "exact.h"

//[[Rcpp::export]]
double get_variance_exact(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const double &tau, const arma::vec &start, const arma::vec &target, 
const arma::vec &lower, const arma::vec &upper, 
const arma::vec &lower_limit, const arma::vec &upper_limit){
  arma::mat Q = get_coffin_matrix(model_name,lower,upper,thetas);
  arma::mat Qt = Q.t();  
  arma::mat large_Q = get_coffin_matrix(model_name,lower_limit,upper_limit,thetas);
  arma::mat large_Qt = large_Q.t();
  arma::mat v(Qt.n_rows,1);
  arma::mat w(1,Q.n_rows);
  arma::mat x(large_Qt.n_rows,1);
  arma::mat y(1,large_Q.n_rows);
  int obs_index = state_to_index(target,lower,upper);
  int start_index = state_to_index(start,lower,upper);
  v[obs_index] = 1;
  arma::mat delta = vT_exp_Q(v,Rcpp::wrap(Qt*tau),1e-20,false,true,false);  // actually Q

  start_index = state_to_index(start,lower_limit,upper_limit);
  obs_index = state_to_index(target,lower_limit,upper_limit);
  y[start_index] = 1;
  x[obs_index] = 1;
  arma::mat a = v_exp_Q(y,Rcpp::wrap(large_Q*(tout-tau)),1e-20,false,true,false); 
  arma::mat b = vT_exp_Q(x,Rcpp::wrap(large_Qt*tau),1e-20,false,true,false); 
  double p = (v_exp_Q(y,Rcpp::wrap(large_Q*tout),1e-20,false,true,false))[obs_index];

  double VRF = p*(1-p);
  arma::vec a_(Q.n_rows-1);
  arma::vec b_(Q.n_rows-1);
  int counter = 0;
  for(int i=0;i<large_Q.n_rows-1;++i){
    if(!(inBox(index_to_state(i,lower_limit,upper_limit),lower,upper))){
      continue;
    }
    a_[counter] = a[i];
    b_[counter] = b[i];
    //std::cout<<"counter:"<<index_to_state(counter,lower,upper)<<" index:"<<index_to_state(i,lower_limit,upper_limit)<<std::endl;
    //std::cout<<"a:"<<a[i]<<", b:"<<b[i]<<", delta"<<b[i]-delta[counter]<<std::endl;
    // percentage_variance_reduction -= (a[i]*b[i]*(1-b[i])-a[i]*(b[i]-delta[counter])*(1-(b[i]-delta[counter]))); 
    VRF -= (a[i]*b[i]*(1-b[i])-a[i]*(b[i]-delta[counter])*(1-(b[i]-delta[counter]))); 
    counter+=1;
  }
  VRF = p*(1-p)/VRF;
  return(VRF);
}
//[[Rcpp::export]]
Rcpp::DataFrame get_variance_exact2(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const double &tau, const arma::vec &start, const arma::vec &target, 
const arma::vec &lower, const arma::vec &upper, 
const arma::vec &lower_limit, const arma::vec &upper_limit){
  arma::mat Q = get_coffin_matrix(model_name,lower,upper,thetas);
  arma::mat Qt = Q.t();  
  arma::mat large_Q = get_coffin_matrix(model_name,lower_limit,upper_limit,thetas);
  arma::mat large_Qt = large_Q.t();
  arma::mat v(Qt.n_rows,1);
  arma::mat w(1,Q.n_rows);
  arma::mat x(large_Qt.n_rows,1);
  arma::mat y(1,large_Q.n_rows);
  int obs_index = state_to_index(target,lower,upper);
  int start_index = state_to_index(start,lower,upper);
  v[obs_index] = 1;
  arma::mat delta = vT_exp_Q(v,Rcpp::wrap(Qt*tau),1e-20,false,true,false);  // actually Q

  start_index = state_to_index(start,lower_limit,upper_limit);
  obs_index = state_to_index(target,lower_limit,upper_limit);
  y[start_index] = 1;
  x[obs_index] = 1;
  arma::mat a = v_exp_Q(y,Rcpp::wrap(large_Q*(tout-tau)),1e-20,false,true,false); 
  arma::mat b = vT_exp_Q(x,Rcpp::wrap(large_Qt*tau),1e-20,false,true,false); 
  double p = (v_exp_Q(y,Rcpp::wrap(large_Q*tout),1e-20,false,true,false))[obs_index];

  double percentage_variance_reduction = p*(1-p);
  arma::vec a2(Q.n_rows-1);
  arma::vec b2(Q.n_rows-1);
  arma::vec delta2(Q.n_rows-1);
  int counter = 0;
  for(int i=0;i<large_Q.n_rows-1;++i){
    if(!(inBox(index_to_state(i,lower_limit,upper_limit),lower,upper))){
      continue;
    }
    a2[counter] = a[i];
    b2[counter] = b[i];
    //std::cout<<"counter:"<<index_to_state(counter,lower,upper)<<" index:"<<index_to_state(i,lower_limit,upper_limit)<<std::endl;
    // std::cout<<index_to_state(counter,lower,upper)<<", a:"<<a[i]<<", b:"<<b[i]<<", delta"<<b[i]-delta[counter]<<std::endl;
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
double get_variance_exact3(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const double &tau, const arma::vec &start, const arma::vec &target, 
const arma::vec &lower, const arma::vec &upper, 
const arma::vec &lower_limit, const arma::vec &upper_limit){
  arma::mat Q = get_coffin_matrix(model_name,lower,upper,thetas);
  arma::mat Qt = Q.t();  
  arma::mat large_Q = get_coffin_matrix(model_name,lower_limit,upper_limit,thetas);
  arma::mat large_Qt = large_Q.t();
  arma::mat v(Qt.n_rows,1);
  arma::mat w(1,Q.n_rows);
  arma::mat x(large_Qt.n_rows,1);
  arma::mat y(1,large_Q.n_rows);
  int obs_index = state_to_index(target,lower,upper);
  int start_index = state_to_index(start,lower,upper);
  v[obs_index] = 1;
  arma::mat delta = arma::expmat(Q*tau);  // actually Q
  start_index = state_to_index(start,lower_limit,upper_limit);
  obs_index = state_to_index(target,lower_limit,upper_limit);
  y[start_index] = 1;
  x[obs_index] = 1;
  arma::mat a = v_exp_Q(y,Rcpp::wrap(large_Q*(tout-tau)),1e-20,false,true,false); 
  arma::mat b = arma::expmat(large_Q*tau); 
  arma::mat p_ = v_exp_Q(y,Rcpp::wrap(large_Q*tout),1e-20,false,true,false);
  double p = p_[obs_index];
  double var_fx = 0;
  arma::vec f(large_Q.n_rows-1);
  arma::vec f2(Q.n_rows-1);
  for(int i=0; i<large_Q.n_rows-1;++i){
    f(i) = get_f(index_to_state(i,lower_limit,upper_limit),target);
  }
  for(int i=0; i<Q.n_rows-1;++i){
    f2(i) = get_f(index_to_state(i,lower,upper),target);
  }
  arma::vec a2(Q.n_rows-1);
  arma::mat b2(Q.n_rows-1,Q.n_rows-1);
  arma::mat delta2(Q.n_rows-1,Q.n_rows-1);
  int i_counter = 0;
  int j_counter = 0;
  
  for(int i=0;i<large_Q.n_rows-1;++i){
    if(!(inBox(index_to_state(i,lower_limit,upper_limit),lower,upper))){
      continue;
    }
    a2[i_counter] = a[i];
    for(int j=0;j<large_Q.n_rows-1;++j){
      if(!(inBox(index_to_state(j,lower_limit,upper_limit),lower,upper))){
        continue;
      }
      b2(i_counter,j_counter) = b(i,j);
      delta2(i_counter,j_counter) = b(i,j)-delta(i_counter,j_counter);
      if(delta2(i_counter,j_counter)<0){
        delta2(i_counter,j_counter) = 0;
      }
      j_counter+=1;
    }
    j_counter=0;
    //std::cout<<index_to_state(counter,lower,upper)<<", "<<delta2[counter]<<std::endl;
    i_counter+=1;
  }
  double component1 = 0;
  double component2 = 0;
  for(int i=0;i<large_Q.n_rows-1;++i){
    component1 += f(i)*f(i)*p_(i);
    component2 += f(i)*p_(i);
  }
  var_fx = component1 - component2*component2;
  // std::cout<<"Var_FX: "<<var_fx<<std::endl;
  double percentage_variance_reduction = var_fx;
  double component3 = 0; 
  double component4 = 0;
  double component5 = 0;
  for(int i=0;i<a2.size();++i){
    for(int j=0;j<a2.size();++j){
      component3 += f2(j)*f2(j)*(b2(i,j)-delta2(i,j));
      component4 += f2(j)*(b2(i,j)-delta2(i,j));
      component5 += f2(j)*(b2(i,j)+delta2(i,j));
    }
    percentage_variance_reduction -= a2(i)*(component3 - component4*component5);
    component3 = 0; 
    component4 = 0;
    component5 = 0;
  }
  percentage_variance_reduction = var_fx/percentage_variance_reduction;
  Rcpp::List df = Rcpp::List::create(
    Rcpp::Named("a")=a2,
    Rcpp::Named("b")=b2,
    Rcpp::Named("delta")=delta2,
    Rcpp::Named("p")=p,
    Rcpp::Named("pvr")=percentage_variance_reduction
  );
  // return(df);
  return(percentage_variance_reduction);
}


//[[Rcpp::export]]
arma::mat get_contour_exact(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const arma::vec &taus, const arma::vec &start, const arma::vec &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list, 
const arma::vec &lower_limit, const arma::vec &upper_limit){
  int number_widths = lower_list.size();
  int number_tau = taus.size();
  arma::vec lower;
  arma::vec upper;
  double tau;
  arma::mat z(number_tau,number_widths);  
  for(int i=0;i<number_widths;++i){
    for(int j=0;j<number_tau;++j){
      lower = Rcpp::as<arma::vec>(lower_list[i]);
      upper = Rcpp::as<arma::vec>(upper_list[i]);
      tau = taus[j];
      z(j,i) = get_variance_exact(model_name, thetas, tout, tau, start, target, lower, upper,lower_limit,upper_limit);
    }
  }
  return(z);
}

//[[Rcpp::export]]
arma::mat get_contour_exact3(const std::string &model_name, const arma::vec &thetas, 
const double &tout, const arma::vec &taus, const arma::vec &start, const arma::vec &target, 
const Rcpp::List &lower_list, const Rcpp::List &upper_list, 
const arma::vec &lower_limit, const arma::vec &upper_limit){
  int number_widths = lower_list.size();
  int number_tau = taus.size();
  arma::vec lower;
  arma::vec upper;
  double tau;
  arma::mat z(number_tau,number_widths);  
  for(int i=0;i<number_widths;++i){
    for(int j=0;j<number_tau;++j){
      lower = Rcpp::as<arma::vec>(lower_list[i]);
      upper = Rcpp::as<arma::vec>(upper_list[i]);
      tau = taus[j];
      z(j,i) = get_variance_exact3(model_name, thetas, tout, tau, start, target, lower, upper,lower_limit,upper_limit);
    }
  }
  return(z);
}

