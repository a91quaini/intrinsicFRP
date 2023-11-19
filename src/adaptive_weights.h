// Author: Alberto Quaini

#ifndef ADAPTIVE_WEIGHTS_H
#define ADAPTIVE_WEIGHTS_H

#include <RcppArmadillo.h>

/////////////////////////////////////////////////////////////
////////  Adaptive Weights //////////////////////////////////
/////////////////////////////////////////////////////////////

// Function for internal use
// Computes adaptive weights from data on factors and test asset
// excess returns. The possible adaptive weights are based on: the first-step
// intrinsic factor risk premia estimates with choice `a`; the matrix of
// factors regression coefficients on test asset excess returns with choice
// `b`; the correlation matrix between factors and returns with choice `c`;
// the unit vector with any other character.
arma::vec AdaptiveWeightsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const char type = 'c'
);

// function for internal use
// Computes n_factors-dimensional vector of adaptive weights based on a
// (n_factors x n_returns)-dimensional matrix, e.g., the correlation matrix
// between factors and returns, the regression coefficients of returns on
// factors (transposed).
arma::vec AdaptiveWeightsFromMatrixCpp(const arma::mat& matrix);

// function for internal use
// Computes n_factors-dimensional vector of adaptive weights based on a
// (n_factors)-dimensional vector, e.g., a first-step intrinsic risk prmeia
// estimator.
arma::vec AdaptiveWeightsFromVectorCpp(const arma::vec& vector);

#endif
