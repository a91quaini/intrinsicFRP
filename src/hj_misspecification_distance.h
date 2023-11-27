// Author: Alberto Quaini

#ifndef HJ_MISSPECIFICATION_DISTANCE_H
#define HJ_MISSPECIFICATION_DISTANCE_H

#include <RcppArmadillo.h>

// HJ asset pricing model misspecification distance.
//
// Computes the Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>
// squared model misspecification distance:
// `square_distance = min_{d} (E[R] - Cov[R,F] * d)' * V[R]^{-1} * (E[R] - Cov[R,F] * d)`,
// where `R` denotes test asset excess returns and `F` risk factors,
// and computes the associated confidence interval.
// This model misspecification distance is a modification of the prominent
// Hansen-Jagannathan (1997) <doi:10.1111/j.1540-6261.1997.tb04813.x>
// distance, adapted to the use of excess returns for the test asset, and a
// SDF that is a linear function of demeaned factors.
// Clearly, computation of the confidence interval is obtained by means of an
// asymptotic analysis under potentially misspecified models, i.e.,
// without assuming correct model specification.
// Details can be found in Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>.
//
// @param returns A `n_observations x n_returns` matrix of test asset excess returns.
// @param factors A `n_observations x n_factors` matrix of risk factors.
// @param ci_coverage A number indicating the confidence interval coverage
// probability. Default is `0.95`.
// @param hac_prewhite A boolean indicating if the series needs pre-whitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// robust covariance (HAC) estimation. Default is `false`.
//
// @return @return A list containing the squared misspecification-robust HJ
// distance in `squared_distance`, and the lower and upper confidence bounds
// in `lower_bound` and `upper_bound`, respectively.
//
// [[Rcpp::export]]
Rcpp::List HJMisspecificationDistanceCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double ci_coverage = .95,
  const bool hac_prewhite = false
);

// Function for internal use
// Computes the Hansen-Jagannatan model misspecification test statistic and p-value
// of an asset pricing model for given test asset excess returns and risk
// factors.
Rcpp::List HJMisspecificationDistanceCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double ci_coverage = .95,
  const bool hac_prewhite = false
);

#endif
