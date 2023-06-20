// Author: Alberto Quaini

#ifndef ADAPTIVE_FRP_H
#define ADAPTIVE_FRP_H

#include "adaptive_ifrp.h"
#include "tuning_adaptive_frp.h"

///////////////////////////////////////////////////
////////  Adaptive Factor Risk Premia /////////////
///////////////////////////////////////////////////

//' Compute optimal adaptive factor risk premia with tuning based on the ratio
//' between model misspecification and mis-identification
//'
//' @name OptimalAdaptiveFRPCpp
//' @description Computes optimal adaptive factor risk premia based
//' on moments extracted from factors and test asset excess returns and adaptive
//' weights over various penalty parameter values. Tuning is performed via
//' minimizing the ratio between model misspecification and mis-identification.
//' Adaptive weights can be based on the
//' correlation between factors and returns, on the regression coefficients of
//' returns on factors or on the first-step intrinsic risk premia estimator.
//'
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of factors.
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and test asset excess returns.
//' @param variance_returns `n_returns x n_returns`-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//' @param penalty_parameters `n_parameters`-dimensional vector of penalty
//' parameter values from smallest to largest.
//' @param weighting_type character specifying the type of adaptive weights:
//' based on the correlation between factors and returns `'c'`; based on the
//' regression coefficients of returns on factors `'b'`; based on the first-step
//' intrinsic risk premia estimator `'a'`; otherwise a vector of ones (any other
//' character). Default is `'c'`.
//'
//' @return a list containing the `n_factors`-dimensional vector of adaptive
//' intrinsic factor risk premia in `"risk_premia"`, and the optimal penalty
//' parameter value in `"penalty_parameter"`.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OptimalAdaptiveFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c'
);

//' Compute the HAC standard errors of the nonzero adaptive factor risk premia
//'
//' @name StandardErrorsAdaptiveFRPCpp
//' @description Computes the HAC standard errors of adaptive factor
//' risk premia based on moments extracted from factors and test asset excess
//' returns. It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param afrp n_factors-dimensional vector of factor risk
//' premia, i.e., the result of function `OptimalAdaptiveFRPCpp`.
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of risk
//' factors.
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and test asset excess returns.
//' @param variance_returns `n_returns x n_returns`-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//' @param mean_factors `n_factors`-dimensional mean vector of risk factors.
//'
//' @noRd
//'
//' @return `n_factors`-dimensional vector of standard errors for adaptive
//' intrinsic factor risk premia.
//'
// [[Rcpp::export]]
arma::vec StandardErrorsAdaptiveFRPCpp(
  const arma::vec& afrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

#endif
