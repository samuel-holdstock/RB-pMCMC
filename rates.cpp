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
