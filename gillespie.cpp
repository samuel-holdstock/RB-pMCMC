#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "gillespie.h"


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
  arma::vec xt_data;
  bool set_xttau = false;
  int counter = 0;
  bool inS=true;
  for(int i=0; i<n_spec; ++i){
    if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
      inS = FALSE;
      break;
    }
  }

  while (tnext<tout) {
    counter+=1;
    for(int i=0; i<n_spec; ++i){
      if(counter>2400){
        // std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
        Rcpp::List results = Rcpp::List::create(Rcpp::Named("xttau_data")=xcurr, Rcpp::Named("inS")=inS, Rcpp::Named("xt_data")=xcurr);
        return(results);
      }
    }
    if (tcurr <= tout-tau && tout-tau < tnext){
      xttau_data = xcurr; // make sure not ref
      set_xttau  = true;
      for(int i=0; i<n_spec; ++i){
        if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
          inS = FALSE;
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
    if(tout-tau <= tcurr){
      for(int i=0; i<n_spec; ++i){
          if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
              inS = FALSE;
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
          inS = FALSE;
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
  bool inS=TRUE;
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
              inS = FALSE;
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
              inS = FALSE;
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
          inS = FALSE;
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
  arma::vec xt_data;
  int counter = 0;
  while (tnext<tout) {
    tcurr=tnext;
    counter+=1;
    for(int i=0; i<n_spec; ++i){
      if(counter>2400){
        // std::cout<<"WARNING: GILLESPIE DIVERGE."<<std::endl;
        Rcpp::List results = Rcpp::List::create(Rcpp::Named("xt_data")=xt_data);
        return(results);
      }
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
