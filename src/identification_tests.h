// Author: Alberto Quaini

#ifndef IDENTIFICATION_TESTS_H
#define IDENTIFICATION_TESTS_H

#include <RcppArmadillo.h>

//' Compute the Chen Fang 2019 beta rank test
//'
//' @name ChenFang2019BetaRankTestCpp
//' @description Computes the Chen fang 2019 rank statistic and p-value of the
//' null that the matrix of regression loadings of test asset excess returns on
//' risk factors has reduced rank. If `target_level_kp2006_rank_test > 0`,
//' it uses the iterative Kleibergen Paap 2006 rank test to estimate the initial
//' rank, with `level = target_level_kp2006_rank_test / n_factors`. If
//' `target_level_kp2006_rank_test <= 0`, the initial rank estimator is taken to
//' be the number of singular values above `n_observations^(-1/4)`. It assumes
//' n_factors < n_returns.
//'
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of risk
//' factors.
//' @param n_bootstrap numeric integer indicating the number of bootstrap
//' samples used to compute the Chen fang 2019 test.
//' @param target_level_kp2006_rank_test numeric level of the Kleibergen Paap
//' 2006 rank test. If it is strictly grater than zero, then the iterative
//' Kleibergen Paap 2006 rank test at `level = target_level_kp2006_rank_test /
//' n_factors` is used to compute an initial estimator of the rank of the factor
//' loadings in the Chen Fang 2019 rank test. Otherwise, the initial rank
//' estimator is taken to be the number of singular values above
//' `n_observations^(-1/4)`. Default is `0.05` (as correction for multiple
//' testing).
//'
//' @return a list containing the Chen fang 2019 rank statistic
//' and the corresponding p-value.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List ChenFang2019BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const unsigned int n_bootstrap = 500,
  const double target_level_kp2006_rank_test = 0.05
);

// For internal use
// Computes the iterative Kleibergen Paap 2006 rank statistics and p-values
// of the test that the matrix of regression loadings of test asset excess
// returns on risk factors has rank q = 0, ..., n_factors - 1.
// It also returns an estimate of the rank as the first value q with associated
// p-value below a given `level = target_level / n_factors`.
// It assumes n_factors < n_returns.
Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level = .05
);

// For internal use
// Computes the Kleibergen Paap 2006 rank statistic and p-value
// of the test that the matrix of regression loadings of test asset excess
// returns on risk factors has rank `q`.
arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::mat& theta_vec,
  const arma::mat& asy_var_theta_vec,
  const arma::mat& U_theta,
  const arma::mat& V_theta,
  const arma::mat& Aq_perp,
  const arma::mat& Bq_perp,
  const unsigned int q
);

// For internal use
// Computes the scaled (n_factors x n_returns) matrix of factor loadings,
// useful for the Kleibergen Paap 2006 or the Chen Fang 2019 rank tests.
// The resulting invariant matrix is proportional to the matrix of t-statistics
// of the least squares regression estimator.
arma::mat ScaledFactorLoadingsCpp(
  const arma::mat& returns,
  const arma::mat& factors
);

#endif
