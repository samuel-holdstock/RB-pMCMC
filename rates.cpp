#include "rates.h"

void RateDictionary::add_rate_function(std::string str, std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> func){
    ratesDictionary[str]=func;
}
std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> RateDictionary::get_rate_function(std::string str){
    return(ratesDictionary[str]);
}

void initialize_rates(RateDictionary &rate_dictionary){
    rate_dictionary.add_rate_function("BDI",BDI_get_rates);
    rate_dictionary.add_rate_function("SIR",SIR_get_rates);
    rate_dictionary.add_rate_function("SEIR",SEIR_get_rates);
    rate_dictionary.add_rate_function("SIR_N",SIR_N_get_rates);
    rate_dictionary.add_rate_function("SEIR_N",SEIR_N_get_rates);    
}

ModelDictionary::ModelDictionary(){
  auto bdi_ptr = std::make_shared<BDI>();
  add_model(bdi_ptr); 
}
void ModelDictionary::add_model(std::shared_ptr<Model> model_ptr){
  model_dictionary[model_ptr->name] = model_ptr;
}
std::shared_ptr<Model> ModelDictionary::get_model(std::string name){
  std::cout<<model_dictionary[name]->name_backup;
  return(model_dictionary[name]);
}


//[[Rcpp::export]]
Rcpp::NumericVector BDI_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) {
  double P=x[0]; // Define species
  double lambda=theta[0], mu=theta[1], gamma=theta[2]; // Define parameters
  Rcpp::NumericVector r(3); // Number of reactions 

  r[0]=lambda*P; // Calculate rates
  r[1]=mu*P;
  r[2]=gamma;

  return r; // Return vector of rates
}

//[[Rcpp::export]]
Rcpp::NumericVector SIR_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) {
  double S=x[0], I=x[1]; // Define species
  //double N=theta[0], beta=theta[1], gamma=theta[2]; // Define parameters
  double beta=theta[0], gamma=theta[1]; // Define parameters
  Rcpp::NumericVector r(2); // Number of reactions

  //r[0]=beta*S*I/N; // Calculate rates
  r[0]=beta*S*I; // Calculate rates
  r[1]=gamma*I;
 
  return r; // Return vector of rates
}

//[[Rcpp::export]]
Rcpp::NumericVector SIR_N_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) {
  double S=x[0], I=x[1]; // Define species
  double N=theta[0], beta=theta[1], gamma=theta[2]; // Define parameters
  Rcpp::NumericVector r(2); // Number of reactions

  r[0]=beta*S*I/N; // Calculate rates
  r[1]=gamma*I;
 
  return r; // Return vector of rates
}

//[[Rcpp::export]]
Rcpp::NumericVector SEIR_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) {
  double S=x[0], E=x[1], I=x[2]; // Define species
  //double N=theta[0], beta=theta[1], a=theta[2], gamma=theta[3]; // Define parameters
  double beta=theta[0], a=theta[1], gamma=theta[2]; // Define parameters
  Rcpp::NumericVector r(3); // Number of reactions 

  //r[0]=beta*S*I/N; // Calculate rates
  r[0]=beta*S*I; // Calculate rates
  r[1]=a*E;
  r[2]=gamma*I;
 
  return r; // Return vector of rates
}

//[[Rcpp::export]]
Rcpp::NumericVector SEIR_N_get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta) {
  double S=x[0], E=x[1], I=x[2]; // Define species
  double N=theta[0], beta=theta[1], a=theta[2], gamma=theta[3]; // Define parameters
  Rcpp::NumericVector r(3); // Number of reactions 

  r[0]=beta*S*I/N; // Calculate rates
  r[1]=a*E;
  r[2]=gamma*I;
 
  return r; // Return vector of rates
}

Rcpp::NumericMatrix Model::populate_matrix(Rcpp::NumericVector values){
  Rcpp::NumericMatrix S(num_states,num_reactions);
  for(int i=0;i<values.length();++i){
    S(i/num_reactions,i%num_reactions) = values[i];
  }
  return(S);
}

Rcpp::NumericVector BDI::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double P=x[0]; // Define species
  double lambda=theta[0], mu=theta[1], gamma=theta[2]; // Define parameters
  Rcpp::NumericVector r(3); // Number of reactions 
  r[0]=lambda*P; // Calculate rates
  r[1]=mu*P;
  r[2]=gamma;
  return r; // Return vector of rates
}
BDI::BDI(){
  name = "BDI";
  name_backup = "BDI";
  num_reactions = 3;
  num_states = 1;
  num_params = 3;
  values = {1,-1,1};
  S = populate_matrix(values);
}

Rcpp::NumericVector SIR::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double S=x[0], I=x[1]; // Define species
  //double N=theta[0], beta=theta[1], gamma=theta[2]; // Define parameters
  double beta=theta[0], gamma=theta[1]; // Define parameters
  Rcpp::NumericVector r(2); // Number of reactions

  //r[0]=beta*S*I/N; // Calculate rates
  r[0]=beta*S*I; // Calculate rates
  r[1]=gamma*I;
 
  return r; // Return vector of rates
}
void SIR::initialize(){
  name = "SIR";
  num_reactions = 2;
  num_states = 2;
  num_params = 2;
  values = {-1,0,
            1,-1};
  S = populate_matrix(values);
}

Rcpp::NumericVector SIR_N::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double S=x[0], I=x[1]; // Define species
  double N=theta[0], beta=theta[1], gamma=theta[2]; // Define parameters
  Rcpp::NumericVector r(2); // Number of reactions

  r[0]=beta*S*I/N; // Calculate rates
  r[1]=gamma*I;
 
  return r; // Return vector of rates
}
void SIR_N::initialize(){
  name = "SIR_N";
  num_reactions = 2;
  num_states = 2;
  num_params = 3;
  values = {-1,0,
            1,-1};
  S = populate_matrix(values);
}

Rcpp::NumericVector SEIR::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double S=x[0], E=x[1], I=x[2]; // Define species
  //double N=theta[0], beta=theta[1], a=theta[2], gamma=theta[3]; // Define parameters
  double beta=theta[0], a=theta[1], gamma=theta[2]; // Define parameters
  Rcpp::NumericVector r(3); // Number of reactions 

  //r[0]=beta*S*I/N; // Calculate rates
  r[0]=beta*S*I; // Calculate rates
  r[1]=a*E;
  r[2]=gamma*I;
 
  return r; // Return vector of rates
}
void SEIR::initialize(){
  name = "SEIR";
  num_reactions = 3;
  num_states = 3;
  num_params = 3;
  values = {-1,0,0,
            1,-1,0,
            0,1,-1};
  S = populate_matrix(values);
}

Rcpp::NumericVector SEIR_N::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double S=x[0], E=x[1], I=x[2]; // Define species
  double N=theta[0], beta=theta[1], a=theta[2], gamma=theta[3]; // Define parameters
  Rcpp::NumericVector r(3); // Number of reactions 

  r[0]=beta*S*I/N; // Calculate rates
  r[1]=a*E;
  r[2]=gamma*I;
 
  return r; // Return vector of rates
}
void SEIR_N::initialize(){
  name = "SEIR_N";
  num_reactions = 3;
  num_states = 3;
  num_params = 4;
  values = {-1,0,0,
            1,-1,0,
            0,1,-1};
  S = populate_matrix(values);
}

