// Author: Alberto Quaini

#ifndef IFRP_H
#define IFRP_H

#include <RcppArmadillo.h>
#include "hac_standard_errors.h"

////////////////////////////////////////////////
////////////  Intrinsic Factor Risk Premia /////
////////////////////////////////////////////////

//' Compute intrinsic factor risk premia from moments extracted from data
//'
//' @name IFRPCpp
//' @description Computes intrinsic factor risk premia based on moments
//' extracted from factors and test asset excess returns.
//'
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and test asset excess returns.
//' @param variance_returns `n_returns x n_returns`-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//'
//' @return `n_factors`-dimensional vector of intrinsic factor risk premia.
//'
//' @noRd
//'
// [[Rcpp::export]]
arma::vec IFRPCpp(
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

//' Compute the HAC standard errors of intrinsic factor risk premia
//'
//' @name StandardErrorsIFRPCpp
//' @description Computes the HAC standard errors of intrinsic factor risk
//' premia based on moments extracted from factors and test asset excess returns.
//' It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param ifrp n_factors-dimensional vector of intrinsic factor risk premia.
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
//' @param mean_factors n_factors-dimensional mean vector of risk factors.
//'
//' @return n_factors-dimensional vector of HAC standard errors for
//' intrinsic factor risk premia.
//'
//' @noRd
//'
// [[Rcpp::export]]
arma::vec StandardErrorsIFRPCpp(
  const arma::vec& ifrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& mean_factors
);

#endif
