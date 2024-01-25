#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include "gillespie.h"


Rcpp::List gillespie_alg(const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, const Rcpp::NumericMatrix &S, double tout,
                            const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau,
                            std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function) {

  int n_spec = S.nrow(); // Number of species
  int n_react = S.ncol(); // Number of reactions
  double rtot, tcurr=0, tnext=0;
  Rcpp::NumericVector xcurr=clone(x0);
  Rcpp::NumericVector r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  Rcpp::NumericVector xttau_data;
  Rcpp::NumericVector xt_data;
  bool inS=true;
  bool set_xttau = false;

  while (tnext<tout) {
    if (tcurr <= tout-tau && tout-tau < tnext){
      xttau_data = clone(xcurr);
      set_xttau  = true;
    }
    tcurr=tnext;
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S(Rcpp::_,i);
            break;
        }
    }
    for(int i=0; i<n_spec; ++i){
        if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
            inS = false;
            break;
        }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  }
  if(set_xttau == false){
    xttau_data = xcurr;
  }
  xt_data = xcurr;
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("xttau_data")=xttau_data, Rcpp::Named("inS")=inS, Rcpp::Named("xt_data")=xt_data);
  return results;
}

Rcpp::List gillespie_alg_entire(const Rcpp::NumericVector &x0, const Rcpp::NumericVector &theta, const Rcpp::NumericMatrix &S, double tout,
                            const Rcpp::NumericVector &lower, const Rcpp::NumericVector &upper, double tau,
                            std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function) {
  int n_spec = S.nrow(); // Number of species
  int n_react = S.ncol(); // Number of reactions
  double rtot, tcurr=0, tnext=0;
  Rcpp::NumericVector xcurr=clone(x0);
  Rcpp::NumericVector r=rates_function(xcurr,theta);
  rtot=sum(r);
  tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  bool inS=TRUE;
  Rcpp::NumericVector xttau_data(n_spec+1);
  Rcpp::NumericMatrix data(0,n_spec+1);
  bool set_xttau = false;
  while (tnext<tout) {
    data = add_row_time(data,xcurr,tcurr);
    if (tcurr <= tout-tau && tout-tau < tnext){
      xttau_data = clone(xcurr);
      set_xttau = true;
    }
    tcurr=tnext;
    double u = R::runif(0,1);
    double cumtot = 0; 
    for(int i=0; i<n_react; ++i){
        cumtot += r[i];
        if(u*rtot<cumtot){
            xcurr += S(Rcpp::_,i);
            break;
        }
    }
    for(int i=0; i<n_spec; ++i){
        inS = TRUE;
        if (lower[i] > xcurr[i] || xcurr[i] > upper[i]){
            inS = FALSE;
            break;
        }
    }
    r=rates_function(xcurr,theta);
    rtot=sum(r);
    tnext=tcurr-log(R::runif(0,1))/rtot; // add Exp(rtot)
  }
  if(set_xttau == false){
    xttau_data = xcurr;
  }
  Rcpp::List results = Rcpp::List::create(Rcpp::Named("xttau_data")=xttau_data, Rcpp::Named("inS")=inS, Rcpp::Named("data")=data);
  return results;
}

Rcpp::NumericMatrix add_row_time(Rcpp::NumericMatrix data, Rcpp::NumericVector xcurr, double tnext){
  int nrows = data.rows();
  int ncols = data.ncol();
  Rcpp::NumericMatrix newData = Rcpp::NumericMatrix(nrows+1,ncols);
  for(int i=0;i<nrows;++i){
    newData(i,Rcpp::_) = data(i,Rcpp::_);
  }
  newData(nrows,0) = tnext;
  for(int j=1;j<=ncols-1;++j){
    newData(nrows,j) = xcurr[j-1];
  }
  return(newData);
}
