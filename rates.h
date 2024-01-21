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

//[[Rcpp::export]]
Rcpp::NumericVector BDI_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);
//[[Rcpp::export]]
Rcpp::NumericVector SIR_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);
//[[Rcpp::export]]
Rcpp::NumericVector SEIR_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);
//[[Rcpp::export]]
Rcpp::NumericVector SIR_N_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);
//[[Rcpp::export]]
Rcpp::NumericVector SEIR_N_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta);

class Model{
  public:
    virtual Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) = 0;
    virtual Rcpp::NumericMatrix get_S();
    int num_reactions;
    int num_states;
    int num_params;
  protected:
    Rcpp::NumericVector values;
    Rcpp::NumericMatrix populate_matrix(Rcpp::NumericVector values);
};

class BDI: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    Rcpp::NumericMatrix get_S() override;
};



#endif