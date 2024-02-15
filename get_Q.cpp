#include "get_Q.h"

// void calculate_Q(Rcpp::NumericVector upper, Rcpp::NumericVector lower){
//     int num_param = upper.length();
//     int num_points = 1;
//     std::vector<int> prod(num_param);
//     int prev = 1;
//     for(int i=0;i<num_param;++i){
//         num_points*=(upper[i]-lower[i]+1);
//         prod[i]=(upper[i]-lower[i]+1)*prev;
//         prev = prod[i];
//         std::cout<<prod[i]<<std::endl;
//     }
//     for(int i=0;i<num_points;++i){
//         std::vector<int> coordinates(num_param);
//         prev = 1;
//         std::cout<<i<<": ";
//         for(int j=0;j<num_param;++j){
//             coordinates[j] = (i/prev) % prod[j]+lower[j];
//             prev = prod[j];
//             std::cout<<coordinates[j]<<" ";
//         }
//         std::cout<<std::endl;
//     }
// } 

//[[Rcpp::export]]
int state_to_index(
Rcpp::NumericVector state,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper){
    int num_param = state.length();
    int prev = 1;
    int index = 0;
    bool coffin = false;
    int state_diff = 0;
    int box_diff = 0;
    for(int i=num_param-1;0<=i;--i){
        state_diff = state[i]-lower[i];
        box_diff = upper[i]-lower[i]+1;
        index += (state_diff % (prev*box_diff)) * prev;
        prev *= box_diff;
        if(state_diff+1>box_diff||state_diff<0){
            coffin = true;
        }
    }
    if(coffin){
        index = prev;
    }
    return(index);
}

//[[Rcpp::export]]
Rcpp::NumericVector index_to_state(
int index,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper){
    int num_param = lower.length();
    int prev = 1;
    int num_states = 0;
    Rcpp::NumericVector coordinates(num_param);
    int index_copy = index;
    int box_diff = 0;
    for(int i=num_param-1;0<=i;--i){
        box_diff = upper[i]-lower[i]+1;
        coordinates[i] = (index % box_diff)+lower[i];
        index /= box_diff;
        num_states += box_diff;
    }
    if(index_copy>=num_states){
        std::cout<<"WARNING INDEX>NUM_STATES"<<std::endl;
    }
    return(coordinates);    
}

bool inBox(
Rcpp::NumericVector state,
Rcpp::NumericVector lower,
Rcpp::NumericVector upper
){
    int num_param = lower.length();
    for(int i=0;i<num_param;++i){
        if(state[i]<lower[i] || state[i]>upper[i]){
            return(FALSE);
        }
    }
    return(TRUE);
}

Rcpp::NumericMatrix Qmat(
Rcpp::NumericVector lower,
Rcpp::NumericVector upper,
Rcpp::NumericVector theta,
Rcpp::NumericMatrix S,
std::function<Rcpp::NumericVector(const Rcpp::NumericVector &x, const Rcpp::NumericVector &theta)> rates_function
){
    int num_param = lower.length();
    int num_react = S.ncol();
    int total_points = 1;
    for(int i=0;i<num_param;++i){
        total_points *= (upper[i]-lower[i]+1);
    }
    Rcpp::NumericVector state(num_param);
    Rcpp::NumericVector temp_state(num_param);
    Rcpp::NumericMatrix Q(total_points+1,total_points+1);
    Rcpp::NumericVector rate(num_react);
    int temp_index;
    double stay_rate;
    int coffin_index = total_points;
    for(int i=0;i<total_points;i++){
        state = index_to_state(i,lower,upper);
        rate = rates_function(state,theta);
        stay_rate = 0;
        for(int j=0;j<num_react;++j){
            stay_rate += rate[j];
            for(int k=0;k<num_param;++k){
                temp_state[k] = state[k]+S(k,j);
            }
            if(inBox(temp_state,lower,upper)){
                temp_index = state_to_index(temp_state,lower,upper);                    
                Q(i,temp_index) = Q(i,temp_index) + rate[j];
            }
            else{
                Q(i,coffin_index) = Q(i,coffin_index) + rate[j];
            }
        }
        Q(i,i) = -stay_rate;
    }
    return(Q);
}

