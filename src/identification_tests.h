// Author: Alberto Quaini

#ifndef IDENTIFICATION_TESTS_H
#define IDENTIFICATION_TESTS_H

#include <RcppArmadillo.h>

// Compute the Chen-Fang (2019) beta rank test.
//
// Computes the Chen-Fang (2019) <doi:10.3982/QE1139> rank statistic and p-value for
// testing whether the matrix of regression loadings of test asset excess returns
// on risk factors has a reduced rank. The Kleibergen-Paap (2006) <doi:10.1016/j.jeconom.2005.02.011>
// iterative rank test is used to estimate the initial rank if `target_level_kp2006_rank_test` is
// greater than 0, adjusting the level as `target_level_kp2006_rank_test / n_factors`.
// If `target_level_kp2006_rank_test` is less than or equal to 0, the initial rank estimator is
// the count of singular values above `n_observations^(-1/4)`. Assumes that `n_factors < n_returns`.
//
// @param returns An `n_observations x n_returns` matrix of test asset excess returns.
// @param factors An `n_observations x n_factors` matrix of risk factors.
// @param n_bootstrap Integer number of bootstrap samples for the Chen-Fang (2019) test.
// @param target_level_kp2006_rank_test Double indicating the Kleibergen-Paap (2006) rank test level.
// If greater than 0, it specifies the level for the initial rank estimation in the
// Chen-Fang (2019) rank test. Otherwise, uses the count of singular values above
// `n_observations^(-1/4)` as the initial rank estimator. Default is `0.05` for
// multiple testing correction.
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
// Computes the iterative Kleibergen-Paap (2006) <doi:10.1016/j.jeconom.2005.02.011> rank statistics and p-values for
// testing if the matrix of regression loadings of test asset excess returns on risk factors
// has rank `q = 0, ..., n_factors - 1`. It also estimates the rank as the first `q` with a p-value
// below the given level, adjusted as `target_level / n_factors`. Assumes `n_factors < n_returns`.
//
// @param returns An `n_observations x n_returns` matrix of test asset excess returns.
// @param factors An `n_observations x n_factors` matrix of risk factors.
// @param target_level Double specifying the target level of the test used for rank estimation.
//
// @return A list containing the iterative Kleibergen-Paap 2006 beta rank statistics and p-values.
// [[Rcpp::export]]
Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level = .05
);

// Compute the Kleibergen-Paap 2006 beta rank test statistic and p-value.
//
// Computes the Kleibergen-Paap (2006) <doi:10.1016/j.jeconom.2005.02.011> rank statistic and p-value for testing if
// the matrix of regression loadings of test asset excess returns on risk factors
// has a specified rank `q`.
//
// @param theta_vectorised Matrix of regression coefficients.
// @param U Left singular vectors from the SVD of the regression coefficients.
// @param V Right singular vectors from the SVD of the regression coefficients.
// @param W Orthogonal complement of the matrix for rank `q`.
// @param q Unsigned integer specifying the hypothesized rank.
//
// @return A vector containing the Kleibergen-Paap 2006 rank statistic and the p-value.
arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::mat& theta_vectorised,
  const arma::mat& U,
  const arma::mat& V,
  const arma::mat& W,
  const unsigned int q
);

// Compute the scaled matrix of factor loadings.
//
// Computes the scaled (n_factors x n_returns) matrix of scaled factor loadings,
// which is useful for the Kleibergen-Paap (2006)
// <doi:10.1016/j.jeconom.2005.02.011> and the
// Chen-Fang (2019) <doi:10.3982/QE1139> rank tests. This
// matrix is proportional to the matrix of t-statistics for the least squares
// regression estimator.
// That is, it is invariant to invertible transformations of the data that are
// identical over all observations.
//
// @param returns An `n_observations x n_returns` matrix of test asset excess returns.
// @param factors An `n_observations x n_factors` matrix of risk factors.
// @return A matrix representing the scaled factor loadings.
arma::mat ScaledFactorLoadingsCpp(
  const arma::mat& returns,
  const arma::mat& factors
);

#endif
