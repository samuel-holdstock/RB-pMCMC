#ifndef GET_Q_H
#define GET_Q_H

#include "main.h"

void calculate_Q(arma::vec upper, arma::vec lower);

int state_to_index(
arma::vec state,
arma::vec lower,
arma::vec upper);

arma::vec index_to_state(
int index,
arma::vec lower,
arma::vec upper);

arma::mat Qmat(
arma::vec lower,
arma::vec upper,
arma::vec theta,
arma::mat S,
std::function<arma::vec(const arma::vec &x, const arma::vec &theta)> rates_function
);

bool inBox(
arma::vec state,
arma::vec lower,
arma::vec upper
);

int get_noise_index(arma::vec noise_lower,arma::vec noise_upper,arma::vec lower,arma::vec upper,arma::vec coords);

#endif