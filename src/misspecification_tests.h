// Author: Alberto Quaini

#ifndef MISSPECIFICATION_TESTS_H
#define MISSPECIFICATION_TESTS_H

#include <RcppArmadillo.h>

// Compute the HJ model misspecification test based on data.
//
// Computes the Hansen-Jagannathan (1997) <doi:10.1111/j.1540-6261.1997.tb04813.x>
// model misspecification statistic and
// p-value of an asset pricing model from test asset excess returns `R` and
// risk factors `F`.
// The statistic is defined as:
// `HJDISTANCE = min_{d} (E[R] - beta * d)' * Var[R]^{-1} * (E[R] - beta * d)`
// where `beta = Cov[R, F] * Var[F]^{-1}` are the regression coefficients of
// test asset excess returns `R` on risk factors `F`.
// Detailed computations and p-value calculations can be found in
// Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>.
//
// @param returns A `n_observations x n_returns` matrix of test asset excess returns.
// @param factors A `n_observations x n_factors` matrix of risk factors.
// @param n_simulations The number of simulations to estimate the p-value.
// @return @return A list containing the squared standardized HJ test statistic and
// the corresponding p-value.
//
// [[Rcpp::export]]
Rcpp::List HJMisspecificationTestCpp(
  const arma::mat& returns,
  const arma::mat& factors
);

// Function for internal use
// Computes the Hansen-Jagannatan misspecification test statistic and p-value
// of an asset pricing model for given test asset excess returns and risk
// factors.
Rcpp::List HJMisspecificationStatisticAndPvalueCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

#endif
