#ifndef RATES_H
#define RATES_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <memory>

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
    Rcpp::NumericMatrix populate_matrix(Rcpp::NumericVector values);
    std::string name;
    std::string name_backup="Test";
    int num_reactions;
    int num_states;
    int num_params;
    Rcpp::NumericMatrix S;
    Rcpp::NumericVector values;
    virtual ~Model() = default;
};

class BDI: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    BDI();
};
class SIR: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    void initialize();
};
class SIR_N: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    void initialize();
};
class SEIR: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    void initialize();
};
class SEIR_N: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    void initialize();
};

class ModelDictionary{
  public:
    ModelDictionary();
    std::map<std::string,std::shared_ptr<Model>> model_dictionary;
    void add_model(std::shared_ptr<Model> model);
    std::shared_ptr<Model> get_model(std::string name);
};


#endif