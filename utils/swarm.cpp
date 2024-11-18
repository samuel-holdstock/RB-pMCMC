#include "swarm.h"

Swarm::Swarm(){
  
}
void Swarm::add_state(const std::pair<arma::vec,int> &state, const double &probability){
  pairs.push_back(std::make_pair(state,probability));
}
void Swarm::push_state(){
  swarm.push_back(pairs);
  pairs.clear();
}
void Swarm::pop_state(){
  swarm.pop_back();
}
std::pair<arma::vec,int> Swarm::get_parent(int parent_index){
  pop_state();
  std::vector<std::pair<std::pair<arma::vec,int>,double>> pairs = swarm.back();
  std::pair<arma::vec,int> parent = pairs[parent_index].first;
  return(parent);
}
arma::mat Swarm::get_path(){
  double normalise_constant = 0;
  std::vector<double> probabilities; 
  std::vector<std::pair<std::pair<arma::vec,int>,double>> pairs = swarm.back();
  for(const auto &entry : pairs){
    normalise_constant += entry.second;
  }
  for(const auto &entry : pairs){
    if(probabilities.empty()){
      probabilities.push_back(entry.second/normalise_constant);
    }
    else{
      probabilities.push_back(probabilities.back() + entry.second/normalise_constant);
    }
  }
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0,1.0);  
  double rand_num;
  int index;
  rand_num = dis(gen);
  auto it = std::lower_bound(probabilities.begin(),probabilities.end(),rand_num);
  index = std::distance(probabilities.begin(),it); 
  arma::mat path(swarm.size(),pairs[index].first.first.size());
  path.row(0)=pairs[index].first.first.t();
  std::pair<arma::vec,int> parent;
  int parent_index = pairs[index].first.second;
  for(int i=0;i<swarm.size();++i){
    parent = get_parent(parent_index);
    path.row(i+1)=parent.first.t();
    parent_index = parent.second;
  }
  return(path);
}
