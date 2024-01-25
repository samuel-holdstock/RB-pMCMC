#ifndef REXPQ_H
#define REXPQ_H

#include <iostream>

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/distributions/poisson.hpp>
#include <chrono>

using namespace std;

inline double h(double x);
inline double hifunc(double rho, double B,double loge);
inline double lofunc(double rho, double A, double loge);
unsigned int get_mlo(unsigned int mhi, double rho);
unsigned int get_m(double rho, double prec, unsigned int mlo=0);
bool commonChecks(double prec, int vn_rows, int vn_cols, int Qn_rows, int Qn_cols, string sinit, bool ignore_v=false) ;
bool sChecks(const arma::mat v, const arma::sp_mat Q, double prec, string caller, bool ignore_v=false);
bool dChecks(const arma::mat v, const arma::dmat Q, double prec, string caller, bool ignore_v=false);
static arma::mat drivel(void );
int get_s_from_rho_prec(double rho, double eps);
void renormalise_P(arma::dmat &P);
void renormalise_v(arma::dmat &v, const arma::dvec &v0sums);
arma::mat sSS_exp_Q(const arma::sp_mat Q, double prec, bool renorm=true);
arma::mat dSS_exp_Q(const arma::dmat Q, double prec, bool renorm=true);

arma::mat SS_exp_Q(SEXP Q, double prec, bool renorm=true);
arma::mat sSS_v_exp_Q(const arma::mat v, const arma::sp_mat Q, double prec,
		      bool renorm=true);
arma::mat dSS_v_exp_Q(const arma::mat v, const arma::mat Q, double prec,
		      bool renorm=true);

arma::mat SS_v_exp_Q(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool checks=true);
arma::mat sUnif_v_exp_Q(const arma::mat v, const arma::sp_mat Q, double prec, bool renorm=true, bool t2=true);
arma::mat dUnif_v_exp_Q(const arma::mat v, const arma::mat Q, double prec, bool renorm=true, bool t2=true);

arma::mat Unif_v_exp_Q(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true);
bool should_use_SS_sparse(int d, double rho, double f);
bool should_use_SS_dense(int d, double rho);
arma::mat sv_exp_Q(const arma::mat v, const arma::sp_mat Q, double prec, bool renorm=true, bool t2=true);
arma::mat dv_exp_Q(const arma::mat v, const arma::mat Q, double prec, bool renorm=true, bool t2=true);

arma::mat v_exp_Q(const arma::mat v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true);

arma::dvec vT_exp_Q(const arma::dvec &v, SEXP Q, double prec, bool renorm=true, bool t2=true, bool checks=true);

#endif