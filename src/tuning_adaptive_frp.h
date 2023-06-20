// Author: Alberto Quaini

#ifndef TUNING_ADAPTIVE_FRP_H
#define TUNING_ADAPTIVE_FRP_H

#include <RcppArmadillo.h>
#include "identification_tests.h"
#include "misspecification_tests.h"

// Function for internal use
// Computes the Generalized Cross Validation score of the adaptive intrinsic
// factor risk premia for each penalty parameter value.
arma::vec MisspecificationIdentificationScoreAdaptiveFRPCpp(
  const arma::mat& aifrp,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

#endif
