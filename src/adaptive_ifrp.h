// Author: Alberto Quaini

#ifndef ADAPTIVE_IFRP_H
#define ADAPTIVE_IFRP_H

#include <RcppArmadillo.h>
#include "ifrp.h"
#include "adaptive_weights.h"
#include "tuning_adaptive_ifrp.h"
#include "hac_standard_errors.h"

/////////////////////////////////////////////////////////////
////////  Adaptive Intrinsic Factor Risk Premia /////////////
/////////////////////////////////////////////////////////////

//' Compute optimal adaptive intrinsic factor risk premia under generalized
//' cross validation
//'
//' @name OptimalAdaptiveIFRPGCVCpp
//' @description Computes optimal adaptive intrinsic factor risk premia based
//' on moments extracted from factors and test asset excess returns and adaptive
//' weights over various penalty parameter values. Tuning is performed via
//' Generalized Cross Validation (GCV). Adaptive weights can be based on the
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
//' @param gcv_vr_weighting boolean `TRUE` for scaling pricing errors by
//' the inverse variance matrix of asset excess returns; `FALSE` otherwise.
//' Default is `FALSE`.
//' @param gcv_aic_scaling (only relevant for `tuning_type ='g'`)
//' boolean `TRUE` for AIC scaling (`1 / n_observations`); `FALSE` for BIC scaling
//' (`log(n_observations) / n_observations`). Default is `TRUE`.
//' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
//' whose score is not higher than one standard error above the score of the
//' best model; `FALSE` for picking the best model. Default is `FALSE`.
//' @param relaxed boolean `TRUE` for re-fitting the model without shrinkage
//' post selection; `FALSE` otherwise. Default is `FALSE`.
//'
//' @return a list containing the `n_factors`-dimensional vector of adaptive
//' intrinsic factor risk premia in `"risk_premia"`, and the optimal penalty
//' parameter value in `"penalty_parameter"`.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OptimalAdaptiveIFRPGCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const bool gcv_vr_weighting = false,
  const bool gcv_aic_scaling = true,
  const bool one_stddev_rule = false,
  const bool relaxed = false
);

//' Compute optimal adaptive intrinsic factor risk premia under cross validation
//'
//' @name OptimalAdaptiveIFRPCVCpp
//' @description Computes optimal adaptive intrinsic factor risk premia based
//' on moments extracted from factors and test asset excess returns and adaptive
//' weights over various penalty parameter values. Tuning is performed via
//' Cross Validation (CV). Adaptive weights can be based on the correlation
//' between factors and returns, on the regression coefficients of returns on
//' factors or on the first-step intrinsic risk premia estimator.
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
//' @param n_folds integer number of k-fold for cross validation. Default is `5`.
//' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
//' whose score is not higher than one standard error above the score of the
//' best model; `FALSE` for picking the best model. Default is `FALSE`.
//' @param relaxed boolean `TRUE` for re-fitting the model without shrinkage
//' post selection; `FALSE` otherwise. Default is `FALSE`.
//'
//' @return a list containing the n_factors-dimensional vector of adaptive
//' intrinsic factor risk premia in "risk_premia", and the optimal penalty
//' parameter value in "penalty_parameter".
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OptimalAdaptiveIFRPCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const unsigned int n_folds = 5,
  const bool one_stddev_rule = false,
  const bool relaxed = false
);

//' Compute optimal adaptive intrinsic factor risk premia under rolling validation
//'
//' @name OptimalAdaptiveIFRPRVCpp
//' @description Computes optimal adaptive intrinsic factor risk premia based
//' on moments extracted from factors and test asset excess returns and adaptive
//' weights over various penalty parameter values. Tuning is performed via
//' Rolling Validation (RV). Adaptive weights can be based on the correlation
//' between factors and returns, on the regression coefficients of returns on
//' factors or on the first-step intrinsic risk premia estimator.
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
//' @param n_train_observations number of observations in the rolling training
//' set. Default is `120`.
//' @param n_test_observations number of observations in the test set. Default
//' is `12`.
//' @param roll_shift number of observation shift when moving from the rolling
//' window to the next one. Default is `12`.
//' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
//' whose score is not higher than one standard error above the score of the
//' best model; `FALSE` for picking the best model. Default is `FALSE`.
//' @param relaxed boolean `TRUE` for re-fitting the model without shrinkage
//' post selection; `FALSE` otherwise. Default is `FALSE`.
//'
//' @return a list containing the n_factors-dimensional vector of adaptive
//' intrinsic factor risk premia in "risk_premia", and the optimal penalty
//' parameter value in "penalty_parameter".
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OptimalAdaptiveIFRPRVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const unsigned int n_train_observations = 120,
  const unsigned int n_test_observations = 12,
  const unsigned int roll_shift = 12,
  const bool one_stddev_rule = false,
  const bool relaxed = false
);


//' Compute adaptive intrinsic factor risk premia
//'
//' @name AdaptiveIFRPCpp
//' @description Computes adaptive intrinsic factor risk premia based on
//' pre-computed intrinsic factor risk premia and adaptive penalty weights for
//' various penalty parameter values.
//'
//' @param ifrp `n_factors`-dimensional vector of intrinsic factor risk premia.
//' Usually, it is the output of function `IFRPCpp`
//' @param weights `n_factors`-dimensional vector of weights for the penalty
//' term. Usually, it is the output of function `AdaptiveWeightsCpp`.
//' @param penalty_parameters `n_parameters`-dimensional vector of penalty
//' parameter values from smallest to largest.
//'
//' @return `n_factors x n_parameters`-dimensional matrix of adaptive
//' intrinsic factor risk premia.
//'
//' @noRd
//'
// [[Rcpp::export]]
arma::mat AdaptiveIFRPCpp(
  const arma::vec& ifrp,
  const arma::vec& weights,
  const arma::vec& penalty_parameter
);

// For internal use
// Computes adaptive intrinsic factor risk premia based on pre-computed
// intrinsic factor risk premia and adaptive penalty weights for a specific
// penalty parameter value.
arma::vec AdaptiveIFRPCpp(
  const arma::vec& ifrp,
  const arma::vec& weights,
  const double penalty_parameter
);

//' Compute adaptive intrinsic factor risk premia
//'
//' @name RelaxedAdaptiveIFRPCpp
//' @description Computes the intrinsic factor risk premia of the factors selected
//' by an adaptive intrinsic factor risk premia.
//'
//' @param aifrp `n_factors`-dimensional vector of adaptive intrinsic factor
//' risk premia.
//' Usually, it is the output of function `AdaptiveIFRPCpp`
//' @param covariance_factors_returns `n_factors x n_returns`-dimensional
//' covariance matrix between factors and test asset excess returns.
//' @param variance_returns `n_returns x n_returns`-dimensional covariance
//' matrix of test asset excess returns.
//' @param mean_returns `n_returns`-dimensional mean vector of test asset excess
//' returns.
//'
//' @return `n_selected_factors x n_parameters`-dimensional matrix of intrinsic
//' factor risk premia.
//'
//' @noRd
//'
// [[Rcpp::export]]
arma::vec RelaxedAdaptiveIFRPCpp(
  const arma::vec& aifrp,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

//' Compute the HAC standard errors of the nonzero adaptive intrinsic factor risk
//' premia
//'
//' @name StandardErrorsAdaptiveIFRPCpp
//' @description Computes the HAC standard errors of adaptive intrinsic factor
//' risk premia based on moments extracted from factors and test asset excess
//' returns. It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param aifrp n_factors-dimensional vector of intrinsic factor risk
//' premia. E.g., the result of function `OptimalAdaptiveIFRPGCV` or
//' `OptimalAdaptiveIFRPCV`.
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
arma::vec StandardErrorsAdaptiveIFRPCpp(
  const arma::vec& aifrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& mean_factors
);

#endif
