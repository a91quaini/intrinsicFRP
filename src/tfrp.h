// Author: Alberto Quaini

#ifndef TFRP_H
#define TFRP_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////
////////////  Tradable Factor Risk Premia //////
////////////////////////////////////////////////

// Compute tradable factor risk premia
//
// @name TFRP
// @description Computes tradable factor risk premia from data on factors `F` and
// test asset excess returns `R`:
// `TFRP = Cov[F, R] * Var[R]^{-1} * E[R]`;
// which are by construction the negative covariance of factors `F` with
// the SDF projection on asset returns, i.e., the minimum variance SDF.
// Optionally computes the corresponding heteroskedasticity and autocorrelation
// robust standard errors using the Newey-West (1994) <doi:10.2307/2297912>
// plug-in procedure to select the number of relevant lags, i.e.,
// `n_lags = 4 * (n_observations/100)^(2/9)`.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>.
//
// @param returns A `n_observations x n_returns`-dimensional matrix of test asset
// excess returns.
// @param factors A `n_observations x n_factors`-dimensional matrix of factors.
// @param include_standard_errors A boolean: `TRUE` if you want to compute the
// tradable factor risk premia HAC standard errors; `FALSE` otherwise. Default
// is `FALSE`.
// @param hac_prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// robust covariance (HAC) estimation. Default is `false`.
// @param check_arguments A boolean: `TRUE` for internal check of all function
// arguments; `FALSE` otherwise. Default is `TRUE`.
//
// @return A list containing `n_factors`-dimensional vector of tradable factor
// risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`, then
// it further includes `n_factors`-dimensional vector of tradable factor risk
// premia standard errors in `"standard_errors"`.
//
// [[Rcpp::export]]
Rcpp::List TFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool include_standard_errors,
  const bool hac_prewhite = false
);

// Function for internal use
// Computes tradable factor risk premia from moments extracted from data
// Details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
arma::vec TFRPCpp(
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// Function for internal use
// Computes the HAC standard errors of tradable factor risk premia from moments
// extracted from data using the Newey-West (1994) <doi:10.2307/2297912>
// plug-in procedure to select the number of relevant lags, i.e.,
// `n_lags = 4 * (n_observations/100)^(2/9)`.
// Details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
arma::vec StandardErrorsTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite = false
);

#endif
