#ifndef IDENTIFICATION_TESTS_H
#define IDENTIFICATION_TESTS_H

#include <RcppArmadillo.h>

// For internal use
// Computes the minimum singular value of a matrix
double MinSingularValue(const arma::mat matrix);

// For internal use
// Computes the Chen fang 2019 rank statistic and p-value of the
// (n_returns x n_factors)-dimensional matrix of regression loadings of
// test asset excess returns on risk factors.
arma::vec2 BetaRankChenFang2019StatisticAndPvalue(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const int sv_threshold_type = 0,
  const unsigned int n_bootstrap = 1000
);

#endif
