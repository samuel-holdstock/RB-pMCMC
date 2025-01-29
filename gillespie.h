#ifndef GILLESPIE_H
#define GILLESPIE_H

#include "main.h"

std::tuple<arma::vec,bool,arma::vec> gillespie_alg2(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);
std::tuple<arma::vec,bool,arma::vec> gillespie_alg_debug(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function,
                            DebugRB &debug);

Rcpp::List gillespie_alg(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

Rcpp::List gillespie_alg_entire(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            const arma::vec &lower, const arma::vec &upper, double tau,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

Rcpp::List gillespie_alg_frac(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

std::tuple<arma::vec> gillespie_alg_frac2(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

Rcpp::List gillespie_alg_entire_frac(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function);

std::tuple<arma::vec> gillespie_alg_frac_debug(const arma::vec &x0, const arma::vec &theta, const arma::mat &S, double tout,
                            std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function, DebugFrac &debug);

arma::mat add_row_time(arma::mat data, arma::vec xcurr, double tnext);

#endif