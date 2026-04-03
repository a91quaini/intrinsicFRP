// Author: Alberto Quaini

#ifndef IDENTIFICATION_TESTS_H
#define IDENTIFICATION_TESTS_H

#include <RcppArmadillo.h>

// Compute the Chen-Fang (2019) beta rank test.
//
// This function implements the full-column-rank case of the Chen-Fang (2019)
// procedure for the matrix of regression loadings of test asset excess returns
// on risk factors. The Kleibergen-Paap (2006) iterative exact-rank test is
// used to estimate the initial rank if `target_level_kp2006_rank_test > 0`,
// adjusting the level as `target_level_kp2006_rank_test / n_factors`.
// Otherwise, the initial rank estimator is the count of singular values above
// `n_observations^(-1/4)`. Assumes `n_factors < n_returns`.
//
// @param returns An `n_observations x n_returns` matrix of test asset excess returns.
// @param factors An `n_observations x n_factors` matrix of risk factors.
// @param n_bootstrap Integer number of bootstrap samples for the Chen-Fang (2019) test.
// @param target_level_kp2006_rank_test Double indicating the Kleibergen-Paap (2006)
//   rank test level used for the initial rank estimation. If greater than 0,
//   the sequential KP estimator uses the first `q` with p-value above
//   `target_level_kp2006_rank_test / n_factors`. Otherwise, the initial rank
//   estimator is the count of singular values above `n_observations^(-1/4)`.
// @return A list containing the Chen-Fang 2019 rank statistic and the corresponding p-value.
// [[Rcpp::export]]
Rcpp::List ChenFang2019BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const unsigned int n_bootstrap = 500,
  const double target_level_kp2006_rank_test = 0.05
);

// Compute the iterative Kleibergen-Paap 2006 beta rank test.
//
// Computes the iterative Kleibergen-Paap (2006) rank statistics and p-values
// for testing whether the matrix of regression loadings of test asset excess
// returns on risk factors has rank `q = 0, ..., n_factors - 1`. It estimates
// the rank as the smallest `q` with p-value above the adjusted level
// `target_level / n_factors`. Assumes `n_factors < n_returns`.
//
// @param returns An `n_observations x n_returns` matrix of test asset excess returns.
// @param factors An `n_observations x n_factors` matrix of risk factors.
// @param target_level Double specifying the target level used for rank estimation.
// @return A list containing the iterative Kleibergen-Paap 2006 beta rank statistics and p-values.
// [[Rcpp::export]]
Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level = .05
);

// Compute the Kleibergen-Paap 2006 beta rank test statistic and p-value.
//
// @param theta_vectorised Vectorised matrix of scaled regression coefficients.
// @param U Left singular vectors from the SVD of the scaled regression coefficients.
// @param V Right singular vectors from the SVD of the scaled regression coefficients.
// @param W Covariance matrix of `vec(theta)`.
// @param q Unsigned integer specifying the hypothesised rank.
// @return A vector containing the Kleibergen-Paap 2006 rank statistic and the p-value.
arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::vec& theta_vectorised,
  const arma::mat& U,
  const arma::mat& V,
  const arma::mat& W,
  const unsigned int q
);

// Compute the scaled matrix of factor loadings.
//
// Computes the scaled `(n_factors x n_returns)` matrix of factor loadings,
// useful for the Kleibergen-Paap (2006) and Chen-Fang (2019) rank tests.
// The resulting invariant matrix is proportional to the matrix of t-statistics
// of the least squares regression estimator.
//
// @param returns An `n_observations x n_returns` matrix of test asset excess returns.
// @param factors An `n_observations x n_factors` matrix of risk factors.
// @return A matrix representing the scaled factor loadings.
arma::mat ScaledFactorLoadingsCpp(
  const arma::mat& returns,
  const arma::mat& factors
);

#endif
