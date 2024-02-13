// Author: Alberto Quaini

#ifndef FGX_THREE_PASS_COVARIANCE_H
#define FGX_THREE_PASS_COVARIANCE_H

#include <RcppArmadillo.h>

// Function for internal use
// Compute the covariance matrix of the "new" factors SDF coefficients in the
// three-pass procedure of Feng Giglio and Xiu (2020)
// <doi:https://doi.org/10.1111/jofi.12883>.
//
// [[Rcpp::export]]
arma::mat FGXThreePassCovarianceCpp(
  const arma::mat& returns,
  const arma::mat& control_factors,
  const arma::mat& new_factors,
  const arma::vec& sdf_coefficients,
  const arma::uvec& idx_selected
);

#endif
