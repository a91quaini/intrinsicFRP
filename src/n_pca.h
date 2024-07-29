// Author: Alberto Quaini

#ifndef N_PCA_H
#define N_PCA_H

#include <RcppArmadillo.h>

// Function for internal use.
// Implement the procedure to determine the number of PCA
// that summarize the risk in returns of Giglio Xiu (2021).
unsigned int NPCA_GiglioXiu2021Cpp(
  const arma::vec& evals,
  const unsigned int n_assets,
  const unsigned int n_observations,
  unsigned int n_max
);

Rcpp::List NPCA_AhnHorenstein2013Cpp(
  const arma::vec& evals,
  unsigned int n_max = 0
);

#endif
