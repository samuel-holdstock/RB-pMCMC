#ifndef RATES_H
#define RATES_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

class RateDictionary{
    public:
        std::map<std::string,std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)>>ratesDictionary;
        void add_rate_function(std::string str,std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> func);
        std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> get_rate_function(std::string str);
};

void initialize_rates(RateDictionary &rate_dictionary);

Rcpp::NumericVector BDI_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);
Rcpp::NumericVector SIR_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);
Rcpp::NumericVector SEIR_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);

#endif