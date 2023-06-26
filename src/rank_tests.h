#ifndef RANK_TESTS_H
#define RANK_TESTS_H

#include <RcppArmadillo.h>

// For internal use
// Computes the Kleibergen Paap 2005 rank statistic and p-value that the
// (n_returns x n_factors)-dimensional matrix of regression loadings of
// test asset excess returns on risk factors has rank q < min(n_returns, n_factors).
arma::vec2 BetaRankKleibergenPaap2006StatisticAndPvalueCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const unsigned int q,
  const bool scaling = true
);

// For internal use
// Computes the Chen fang 2019 rank statistic and p-value of the
// (n_returns x n_factors)-dimensional matrix of regression loadings of
// test asset excess returns on risk factors.
arma::vec2 BetaRankChenFang2019StatisticAndPvalueCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const int sv_threshold_type = 0,
  const unsigned int n_bootstrap = 500
);

#endif
