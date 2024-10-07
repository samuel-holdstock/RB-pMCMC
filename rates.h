#ifndef RATES_H
#define RATES_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <memory>

class Model{
  public:
    virtual arma::vec get_rates(const arma::vec &x, const arma::vec &theta) = 0;
    //std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> get_rate_function();
    arma::mat populate_matrix(arma::vec values);
    std::string name;
    int num_reactions;
    int num_states;
    int num_params;
    arma::mat S;
    arma::vec values;
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
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    BDI();
};
class SIR: public Model{
  public:
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    SIR();
};
class SIR_N: public Model{
  public:
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    SIR_N();
};
class SEIR: public Model{
  public:
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    SEIR();
};
class SEIR_N: public Model{
  public:
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    SEIR_N();
};
class BDI_2: public Model{
  public:
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    BDI_2();
};
class LV: public Model{
  public:
    arma::vec get_rates(const arma::vec &x, const arma::vec &theta) override;
    LV();
};

extern ModelDictionary model_dict;

#endif