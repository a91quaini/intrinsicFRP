// Author: Alberto Quaini

#ifndef TFRP_H
#define TFRP_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////
////////////  Tradable Factor Risk Premia //////
////////////////////////////////////////////////

//' Compute tradable factor risk premia from data
//'
//' @name TFRPCpp
//' @description Computes tradable factor risk premia from data on factors `F` and
//' test asset excess returns `R`:
//' `TFRP = Cov[F, R] * Var[R]^{-1} * E[R]`.
//' Optionally computes the corresponding heteroskedasticity and autocorrelation
//' robust standard errors using the Newey-West (1994) <doi:10.2307/2297912>
//' plug-in procedure to select the number of relevant lags, i.e.,
//' `n_lags = 4 * (n_observations/100)^(2/9)`.
//' All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>.
//'
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of factors.
//' @param include_standard_errors boolean `TRUE` if you want to compute the
//' tradable factor risk premia HAC standard errors; `FALSE` otherwise. Default
//' is `FALSE`.
//'
//' @return a list containing `n_factors`-dimensional vector of tradable factor
//' risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`, then
//' it further includes `n_factors`-dimensional vector of tradable factor risk
//' premia standard errors in `"standard_errors"`.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List TFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool include_standard_errors
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
  const arma::vec& mean_returns
);

#endif
