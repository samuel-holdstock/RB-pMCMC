#include "gillespie.h"

std::tuple<arma::vec,bool,arma::vec> gillespie_alg2(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function) {

  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr(n_spec);
  xcurr=x0; // make sure not pass by ref
  arma::vec r(n_react);
  r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  // tnext=tcurr-log(dis(gen))/rtot; // add Exp(rtot)
  arma::vec xttau_data(n_spec);
  arma::vec xt_data(n_spec);
  xt_data=x0;
  bool set_xttau = false;
  int counter = 0;
  bool inS=true;
  double u;
  double cumtot;

  while (tnext<tout) {
    counter+=1;
    if(counter>3000){
      std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
      std::tuple<arma::vec,bool,arma::vec> results = std::make_tuple(xttau_data,inS,xt_data);
      return(results);
    }
    if (tout-tau < tnext && tcurr <= tout-tau){
      xttau_data = xcurr; // make sure not ref
      set_xttau  = true;
      for(int i=0; i<n_spec; ++i){
        if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = false;
          break;
        }
      }
    }
    tcurr=tnext;
    u = R::runif(0,1);
    cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    if(inS && (tout-tau <= tcurr)){
      for(int i=0; i<n_spec; ++i){
          if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
              inS = false;
              break;
          }
      }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
    // tnext=tcurr-log(dis(gen))/rtot; // add Exp(rtot)
  }
  if(set_xttau == false){
    xttau_data = xcurr;
    for(int i=0; i<n_spec; ++i){
      if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = false;
          break;
      }
    }
  }
  xt_data = xcurr;  
  std::tuple<arma::vec,bool,arma::vec> results = std::make_tuple(xttau_data,inS,xt_data);
  return results;
}
std::tuple<arma::vec,bool,arma::vec> gillespie_alg_debug(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function, DebugRB &debug) {
  uint64_t before_init = timeSinceEpochNanosec();
  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr(n_spec);
  xcurr=x0; // make sure not pass by ref
  arma::vec r(n_react);
  r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  // tnext=tcurr-log(dis(gen))/rtot; // add Exp(rtot)
  arma::vec xttau_data(n_spec);
  arma::vec xt_data(n_spec);
  xt_data=x0;
  bool set_xttau = false;
  int counter = 0;
  bool inS=true;
  double u;
  double cumtot;
  uint64_t after_init = timeSinceEpochNanosec();
  uint64_t before_enter_box_check;
  uint64_t after_enter_box_check;
  uint64_t after_reaction;
  uint64_t after_left_box_check;
  while (tnext<tout) {
    counter+=1;
    if(counter>3000){
      std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
      std::tuple<arma::vec,bool,arma::vec> results = std::make_tuple(xttau_data,inS,xt_data);
      return(results);
    }
    before_enter_box_check = timeSinceEpochNanosec();
    if ((tout-tau < tnext) && inS && (tcurr <= tout-tau)){
      xttau_data = xcurr; // make sure not ref
      set_xttau  = true;
      for(int i=0; i<n_spec; ++i){
        if ((lower[i] > xcurr[i]) || (xcurr[i] > upper[i])){
          inS = false;
          break;
        }
      }
    }
    after_enter_box_check = timeSinceEpochNanosec();
    debug.enter_box_check_time += after_enter_box_check-before_enter_box_check;
    tcurr=tnext;
    u = R::runif(0,1);
    cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    after_reaction = timeSinceEpochNanosec();
    debug.reaction_time += after_reaction-after_enter_box_check;
    if((tout-tau <= tcurr) && inS){
      for(int i=0; i<n_spec; ++i){
          if ((lower[i] > xcurr[i]) || (xcurr[i] > upper[i])){
              inS = false;
              break;
          }
      }
    }
    after_left_box_check = timeSinceEpochNanosec();
    debug.left_box_check_time += after_left_box_check-after_reaction;
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
    // tnext=tcurr-log(dis(gen))/rtot; // add Exp(rtot)
  }
  u_int64_t before_no_reaction_check = timeSinceEpochNanosec();
  if(set_xttau == false){
    xttau_data = xcurr;
    for(int i=0; i<n_spec; ++i){
      if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = false;
          break;
      }
    }
  }
  u_int64_t after_no_reaction_check = timeSinceEpochNanosec();
  xt_data = xcurr;  
  std::tuple<arma::vec,bool,arma::vec> results = std::make_tuple(xttau_data,inS,xt_data);
  uint64_t end = timeSinceEpochNanosec();
  debug.init_time += after_init-before_init;
  debug.no_reaction_check_time += after_no_reaction_check-before_no_reaction_check; 
  return results;
}


Rcpp::List gillespie_alg(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function) {

  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr=x0; // make sure not pass by ref
  arma::vec r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  arma::vec xttau_data;
  arma::vec xt_data=x0;
  bool set_xttau = false;
  int counter = 0;
  bool inS=true;
  while (tnext<tout) {
    counter+=1;
    if(counter>3000){
      std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
      Rcpp::List results = Rcpp::List::create(Rcpp::Named("xttau_data")=xcurr, Rcpp::Named("inS")=inS, Rcpp::Named("xt_data")=xcurr);
      return(results);
    }
    if (tcurr <= tout-tau && tout-tau < tnext){
      xttau_data = xcurr; // make sure not ref
      set_xttau  = true;
      for(int i=0; i<n_spec; ++i){
        if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = false;
          break;
        }
      }
    }
    tcurr=tnext;
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    if(inS && (tout-tau <= tcurr)){
      for(int i=0; i<n_spec; ++i){
          if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
              inS = false;
              break;
          }
      }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  }
  if(set_xttau == false){
    xttau_data = xcurr;
    for(int i=0; i<n_spec; ++i){
      if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = false;
          break;
      }
    }
  }
  xt_data = xcurr;
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("xttau_data")=xttau_data, Rcpp::Named("inS")=inS, Rcpp::Named("xt_data")=xt_data);
  return results;
}

Rcpp::List gillespie_alg_entire(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function) {
  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr=x0; // check not ref
  arma::vec r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  bool inS=true;
  arma::vec xttau_data(n_spec+1);
  arma::mat data(0,n_spec+1);
  bool set_xttau = false;

  data = add_row_time(data,xcurr,tcurr);
  while (tnext<tout) {
    if (tcurr <= tout-tau && tout-tau < tnext){
      xttau_data = xcurr; // check not ref
      set_xttau = true;
      for(int i=0; i<n_spec; ++i){
          if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
              inS = false;
              break;
          }
      }
    }
    tcurr=tnext;
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    data = add_row_time(data,xcurr,tcurr);
    if(tout-tau <= tcurr){
      for(int i=0; i<n_spec; ++i){
          if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
              inS = false;
              break;
          }
      }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  }
  if(set_xttau == false){
    xttau_data = xcurr;
    for(int i=0; i<n_spec; ++i){
      if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = false;
          break;
      }
    }
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("xttau_data")=xttau_data, Rcpp::Named("inS")=inS, Rcpp::Named("data")=data, Rcpp::Named("xt_data")=xcurr);
  return results;
}

Rcpp::List gillespie_alg_frac(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function) {
  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr=x0; // clone
  arma::vec r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  arma::vec xt_data=x0;
  int counter = 0;
  while (tnext<tout) {
    tcurr=tnext;
    counter+=1;
    if(counter>3000){
      std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
      Rcpp::List results = Rcpp::List::create(Rcpp::Named("xt_data")=xt_data);
      return(results);
    }
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  xt_data = xcurr;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("xt_data")=xt_data);
  return results;
}
std::tuple<arma::vec> gillespie_alg_frac2(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function) {
  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr=x0; // clone
  arma::vec r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  arma::vec xt_data=x0;
  int counter = 0;
  while (tnext<tout) {
    tcurr=tnext;
    counter+=1;
    if(counter>3000){
      std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
      std::tuple<arma::vec> results = std::make_tuple(xt_data);
      return(results);
    }
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  xt_data = xcurr;
  }
  std::tuple<arma::vec> results = std::make_tuple(xt_data);
  return results;
}
std::tuple<arma::vec> gillespie_alg_frac_debug(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function, DebugFrac &debug){

  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr=x0; // clone
  arma::vec r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  arma::vec xt_data=x0;
  int counter = 0;
  while (tnext<tout) {
    tcurr=tnext;
    counter+=1;
    if(counter>3000){
      std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
      std::tuple<arma::vec> results = std::make_tuple(xt_data);
      return(results);
    }
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  xt_data = xcurr;
  }
  std::tuple<arma::vec> results = std::make_tuple(xt_data);
  return results;
}

Rcpp::List gillespie_alg_entire_frac(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function) {
  int n_spec = S.n_rows; // Number of species
  int n_react = S.n_cols; // Number of reactions
  double rtot, tcurr=0, tnext=0;
  arma::vec xcurr=x0; // clone
  arma::vec r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  arma::mat data(0,n_spec+1);

  data = add_row_time(data,xcurr,tcurr);
  while (tnext<tout) {
    tcurr=tnext;
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S.col(i);
            break;
        }
    }
    data = add_row_time(data,xcurr,tcurr);
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("data")=data, Rcpp::Named("xt_data")=xcurr);
  return results;
}

arma::mat add_row_time(arma::mat data, arma::vec xcurr, double tnext){
  int nrows = data.n_rows;
  int ncols = data.n_cols;
  arma::mat newData = arma::mat(nrows+1,ncols);
  for(int i=0;i<nrows;++i){
    newData.row(i) = data.row(i);
  }
  newData(nrows,0) = tnext;
  for(int j=1;j<=ncols-1;++j){
    newData(nrows,j) = xcurr[j-1];
  }
  return(newData);
}
