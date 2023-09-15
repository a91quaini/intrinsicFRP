// Author: Alberto Quaini

#ifndef MISSPECIFICATION_TESTS_H
#define MISSPECIFICATION_TESTS_H

#include <RcppArmadillo.h>

//' Compute the Hansen-Jagannatan misspecification test
//'
//' @name HJMisspecificationTestCpp
//' @description Computes the Hansen-Jagannatan misspecification statistic and
//' p-value of an asset pricing model from test asset excess returns and
//' risk factors.
//'
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of risk
//' factors.
//'
//' @return a list containing the HJ test statistic and the corresponding
//' p-value.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List HJMisspecificationTestCpp(
  const arma::mat& returns,
  const arma::mat& factors
);

// For internal use
// Computes the Hansen-Jagannatan misspecification test statistic and p-value
// of an asset pricing model for given test asset excess returns and risk
// factors.
Rcpp::List HJMisspecificationStatisticAndPvalueCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// // For internal use
// // Computes the Hansen-Jagannatan distance of an
// // asset pricing model for given test asset excess returns and risk factors.
// double HJDistanceCpp(
//   const arma::mat& beta,
//   const arma::mat& variance_returns,
//   const arma::vec& mean_returns
// );

#endif
