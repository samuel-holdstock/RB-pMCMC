#include "rates.h"

ModelDictionary model_dict;

ModelDictionary::ModelDictionary(){
  auto BDI_ptr = std::make_shared<BDI>();
  auto SIR_ptr = std::make_shared<SIR>();
  auto SIR_N_ptr = std::make_shared<SIR_N>();
  auto SEIR_ptr = std::make_shared<SEIR>();
  auto SEIR_N_ptr = std::make_shared<SEIR_N>();
  auto BDI_2_ptr = std::make_shared<BDI_2>();
  auto LV_ptr = std::make_shared<LV>();
  
  add_model(BDI_ptr);
  add_model(SIR_ptr);
  add_model(SIR_N_ptr);
  add_model(SEIR_ptr);
  add_model(SEIR_N_ptr);
  add_model(BDI_2_ptr);
  add_model(LV_ptr);
}
void ModelDictionary::add_model(std::shared_ptr<Model> model_ptr){
  model_dictionary[model_ptr->name] = model_ptr;
}
std::shared_ptr<Model> ModelDictionary::get_model(std::string name){
  return(model_dictionary[name]);
}

Rcpp::NumericMatrix Model::populate_matrix(Rcpp::NumericVector values){
  Rcpp::NumericMatrix S(num_states,num_reactions);
  for(int i=0;i<values.length();++i){
    S(i/num_reactions,i%num_reactions) = values[i];
  }
  return(S);
}

// std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> Model::get_rate_function(){
//   return[this](const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
//     return(get_rates(x, theta));
//   };
//}

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
SIR::SIR(){
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
SIR_N::SIR_N(){
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
SEIR::SEIR(){
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
SEIR_N::SEIR_N(){
  name = "SEIR_N";
  num_reactions = 3;
  num_states = 3;
  num_params = 4;
  values = {-1,0,0,
            1,-1,0,
            0,1,-1};
  S = populate_matrix(values);
}

Rcpp::NumericVector BDI_2::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double P=x[0]; // Define species
  double gamma=theta[0], epsilon=theta[1]; // Define parameters
  Rcpp::NumericVector r(2); // Number of reactions 
  r[0]=gamma*P; // Calculate rates
  r[1]=epsilon*P;
  return r; // Return vector of rates
}
BDI_2::BDI_2(){
  name = "BDI_2";
  num_reactions = 2;
  num_states = 1;
  num_params = 2;
  values = {1,-1};
  S = populate_matrix(values);
}

Rcpp::NumericVector LV::get_rates(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta){
  double predator=x[0], prey=x[1]; // Define species
  double alpha=theta[0], beta=theta[1], gamma=theta[2]; // Define parameters
  // alpha is the prey reproduction rate
  // beta is the predation rate
  // gamma is the predator death rate
  Rcpp::NumericVector r(3); // Number of reactions 
  r[0]=alpha*prey; // Calculate rates
  r[1]=beta*prey*predator;
  r[2]=gamma*predator;
  return r; // Return vector of rates
}
LV::LV(){
  name = "LV";
  num_reactions = 3;
  num_states = 2;
  num_params = 3;
  values = {0, 1,-1, 
            1,-1, 0};
  S = populate_matrix(values);
}
