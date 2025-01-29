#include "debug.h"

Rcpp::List DebugRB::get_tuple(){
  Rcpp::List tupleList = Rcpp::List::create(
    Rcpp::_["x"] = x,
    Rcpp::_["log_prob"] = log_prob,
    Rcpp::_["inS"] = inS,
    Rcpp::_["hitS"] = hitS,
    Rcpp::_["xt_data"] = xt_data,
    Rcpp::_["xttau_data"] = xttau_data,
    Rcpp::_["presample"] = presample,
    Rcpp::_["weights"] = weights,
    Rcpp::_["particle_weights"] = particle_weights,
    Rcpp::_["exp_time"] = exp_time,
    Rcpp::_["entire_time"] = entire_time,
    Rcpp::_["gillespie_time"] = gillespie_time,
    Rcpp::_["weight_time"] = weight_time,
    Rcpp::_["resample_time"] = resample_time,
    Rcpp::_["hit_box_check_time"] = enter_box_check_time,
    Rcpp::_["total_reaction_time"] = reaction_time,
    Rcpp::_["in_box_check_time"] = left_box_check_time,
    Rcpp::_["no_reaction_check_time"] = no_reaction_check_time,
    Rcpp::_["init_time"] = init_time
  );
  return(tupleList);
}


Rcpp::List DebugFrac::get_tuple(){
  Rcpp::List tupleList = Rcpp::List::create(
    Rcpp::_["x"] = x,
    Rcpp::_["log_prob"] = log_prob,
    Rcpp::_["presample"] = presample,
    Rcpp::_["weights"] = weights,
    Rcpp::_["entire_time"] = entire_time,
    Rcpp::_["gillespie_time"] = gillespie_time,
    Rcpp::_["weight_time"] = weight_time,
    Rcpp::_["resample_time"] = resample_time
  );
  return(tupleList);  
};