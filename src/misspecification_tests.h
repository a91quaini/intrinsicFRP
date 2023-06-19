// Author: Alberto Quaini

#ifndef MISSPECIFICATION_TESTS_H
#define MISSPECIFICATION_TESTS_H

#include <RcppArmadillo.h>
#include "hac_standard_errors.h"

// For internal use
// Computes the Hansen-Jagannatan distance of an
// asset pricing model for given test asset excess returns and risk factors.
double HJDistance(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// For internal use
// Computes the Hansen-Jagannatan misspecification statistic and p-value of an
// asset pricing model for given test asset excess returns and risk factors.
arma::vec2 HJMisspecificationStatisticAndPvalue(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

#endif
