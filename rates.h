#ifndef RATES_H
#define RATES_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <memory>

class Model{
  public:
    virtual Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) = 0;
    //std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> get_rate_function();
    Rcpp::NumericMatrix populate_matrix(Rcpp::NumericVector values);
    std::string name;
    int num_reactions;
    int num_states;
    int num_params;
    Rcpp::NumericMatrix S;
    Rcpp::NumericVector values;
    virtual ~Model() = default;
};

class ModelDictionary{
  public:
    ModelDictionary();
    std::map<std::string,std::shared_ptr<Model>> model_dictionary;
    void add_model(std::shared_ptr<Model> model);
    std::shared_ptr<Model> get_model(std::string name);
};

class BDI: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    BDI();
};
class SIR: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    SIR();
};
class SIR_N: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    SIR_N();
};
class SEIR: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    SEIR();
};
class SEIR_N: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    SEIR_N();
};
class BDI_2: public Model{
  public:
    Rcpp::NumericVector get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) override;
    BDI_2();
};

extern ModelDictionary model_dict;

#endif