// Author: Alberto Quaini

#ifndef FRP_H
#define FRP_H

#include <RcppArmadillo.h>

/////////////////////////////////////////////////
////////////  Factor Risk Premia ////////////////
/////////////////////////////////////////////////

//' Compute factor risk premia
//'
//' @name FRPCpp
//' @description Computes Fama MacBeth (1973) factor risk premia or misspecification-robust factor
//' risk premia of Kan Robotti Shanken (2013) from data on test asset excess
//' returns and risk factors. Optionally computes the corresponding
//' heteroskedasticity and autocorrelation robust standard errors using the
//' Newey-West (1994) plug-in procedure to select the number of relevant lags,
//' i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of risk
//' factors.
//' @param misspecification_robust boolean `TRUE` for the
//' "misspecification-robust" Kar Robotti Shanken GLS approach using the
//' inverse covariance matrix of returns; `FALSE` for standard Fama-Mac-Beth
//' risk premia. Default is `TRUE`.
//' @param include_standard_errors boolean `TRUE` if you want to compute the
//' adaptive intrinsic factor risk premia HAC standard errors; `FALSE` otherwise.
//' Default is `FALSE`.
//'
//' @return a list containing `n_factors`-dimensional vector of factor
//' risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`, then
//' it further includes `n_factors`-dimensional vector of factor risk
//' premia standard errors in `"standard_errors"`.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List FRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust = true,
  const bool include_standard_errors = false
);

//' Compute Fama-MacBeth (1973) factor risk premia from moments extracted from data
//'
//' @name FMFRPCpp
//' @description Computes Fama-MacBeth (1973) factor risk premia based on moments
//' extracted from factors and test asset returns.
//'
//' @param beta `n_returns x n_factors`-dimensional regression coefficient
//' matrix of test asset returns on risk factors: `beta =
//' covariance(returns, factors) * variance(factors)^(-1)`.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//'
//' @return `n_factors`-dimensional vector of factor risk premia.
//'
//' @noRd
arma::vec FMFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns
);

//' Compute the misspecification robust factor risk premia of Kan Robotti
//' Shanken (2013) from moments extracted from data
//'
//' @name KRSFRP
//' @description Computes misspecification-robust factor risk premia of
//' Kan Robotti Shanken (2013) based on moments extracted from factors and test
//' asset excess returns.
//'
//' @param beta `n_returns x n_factors`-dimensional regression coefficient
//' matrix of test asset returns on risk factors: `beta =
//' covariance(returns, factors) * variance(factors)^(-1)`.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//' @param weighting_matrix `n_returns x n_returns`-dimensional weighting
//' matrix. Typically it is the inverse covariance matrix or the
//' second-moment matrix of test asset excess returns, or the asymptotic
//' covariance of the pricing errors.
//'
//' @return `n_factors`-dimensional vector of factor risk premia.
//'
//' @noRd
arma::vec KRSFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns,
  const arma::mat& weighting_matrix
);

//' Compute the standard errors of Fama-MacBeth factor risk premia from moments
//'
//' @name StandardErrorsFRPCpp
//' @description Computes the HAC standard errors of Fama MacBeth (1973) factor risk
//' premia based on moments extracted from factors and test asset excess returns.
//' It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param frp `n_factors`-dimensional vector of factor risk premia.
//' E.g., the result of function `FRPCpp`.
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of risk
//' factors.
//' @param beta `n_returns x n_factors`-dimensional regression coefficient
//' matrix of test asset returns on risk factors: `beta =
//' covariance(returns, factors) * variance(factors)^(-1)`.
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//'
//'
//' @return `n_factors`-dimensional vector of standard errors for the Fama-MacBeth
//' factor risk premia.
//'
//' @noRd
arma::vec StandardErrorsFRPCpp(
  const arma::vec& frp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

//' Compute the standard errors of krs factor risk premia from moments
//'
//' @name StandardErrorsKRSFRPCpp
//' @description Computes the HAC standard errors of misspecification-robust factor
//' risk premia of Kan Robotti Shanken (2013) based on moments extracted from
//' factors and test asset excess returns. It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param krs_frp `n_factors`-dimensional vector of krs factor risk premia.
//' E.g., the result of function `KRSFRPCpp`.
//' @param returns `n_observations x n_returns`-dimensional matrix of test asset
//' excess returns.
//' @param factors `n_observations x n_factors`-dimensional matrix of risk
//' factors.
//' @param beta `n_returns x n_factors`-dimensional regression coefficient
//' matrix of test asset returns on risk factors: `beta =
//' covariance(returns, factors) * variance(factors)^(-1)`.
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and  test asset excess returns.
//' @param variance_returns `n_returns x n_returns`-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//'
//' @return `n_factors`-dimensional vector of HAC standard errors for the krs factor
//' risk premia.
//'
//' @noRd
arma::vec StandardErrorsKRSFRPCpp(
  const arma::vec& krs_frp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

#endif
