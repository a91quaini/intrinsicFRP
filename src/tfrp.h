// Author: Alberto Quaini

#ifndef TFRP_H
#define TFRP_H

#include <RcppArmadillo.h>

////////////////////////////////////////////////
////////////  Tradable Factor Risk Premia //////
////////////////////////////////////////////////

//' Compute tradable factor risk premia and their standard errors
//'
//' @name TFRPCpp
//' @description Computes tradable factor risk premia based on data on
//' risk factors and test asset excess returns.
//'
//' @param returns `n_observations x n_returns`-dimensional matrix of test
//' asset excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of factors.
//' @param include_standard_errors boolean `TRUE` if you want to compute the
//' tradable factor risk premia HAC standard errors; `FALSE` otherwise.
//' Default is `FALSE`.
//'
//' @return a list containing the `n_factors`-dimensional vector of tradable
//' factor risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`,
//' then it further includes `n_factors`-dimensional vector of factor risk
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

//' Compute tradable factor risk premia from moments extracted from data
//'
//' @name TFRPCpp
//' @description Computes tradable factor risk premia based on moments
//' extracted from factors and test asset excess returns.
//'
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and test asset excess returns.
//' @param variance_returns `n_returns x n_returns`-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//'
//' @return `n_factors`-dimensional vector of tradable factor risk premia.
//'
//' @noRd
arma::vec TFRPCpp(
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

//' Compute the HAC standard errors of tradable factor risk premia
//'
//' @name StandardErrorsTFRPCpp
//' @description Computes the HAC standard errors of tradable factor risk
//' premia based on moments extracted from factors and test asset excess returns.
//' It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param ifrp n_factors-dimensional vector of tradable factor risk premia.
//' E.g., the result of function 'IFRPFCpp'.
//' @param returns (n_observations x n_returns)-dimensional matrix of test asset
//' excess returns.
//' @param factors (n_observations x n_factors)-dimensional matrix of risk
//' factors.
//' @param covariance_factors_returns (n_factors x n_returns)-dimensional
//' covariance matrix between factors and  test asset excess returns.
//' @param variance_returns (n_returns x n_returns)-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns n_returns-dimensional mean vector of test asset excess
//' returns.
//'
//' @return n_factors-dimensional vector of HAC standard errors for
//' tradable factor risk premia.
//'
//' @noRd
arma::vec StandardErrorsTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

#endif
