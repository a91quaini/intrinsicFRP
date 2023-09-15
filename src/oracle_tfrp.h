// Author: Alberto Quaini

#ifndef ORACLE_TFRP_H
#define ORACLE_TFRP_H

#include <RcppArmadillo.h>

/////////////////////////////////////////////////////////////
////////  Oracle Tradable Factor Risk Premia /////////////
/////////////////////////////////////////////////////////////

//' Compute optimal oracle tradable factor risk premia under generalized
//' cross validation
//'
//' @name OracleTFRPGCVCpp
//' @description Computes optimal oracle tradable factor risk premia based
//' on moments extracted from factors and test asset excess returns and oracle
//' weights over various penalty parameter values. Tuning is performed via
//' Generalized Cross Validation (GCV). Oracle weights can be based on the
//' correlation between factors and returns, on the regression coefficients of
//' returns on factors or on the first-step tradable risk premia estimator.
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
//' @param weighting_type character specifying the type of oracle weights:
//' based on the correlation between factors and returns `'c'`; based on the
//' regression coefficients of returns on factors `'b'`; based on the first-step
//' tradable risk premia estimator `'a'`; otherwise a vector of ones (any other
//' character). Default is `'c'`.
//' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
//' whose score is not higher than one standard error above the score of the
//' best model; `FALSE` for picking the best model. Default is `FALSE`.
//' @param gcv_scaling_n_assets boolean `TRUE` for sqrt(n_assets) scaling
//' (`sqrt(n_assets) / n_observations`); `FALSE` otherwise (`1 / n_observations`).
//' Default is `FALSE`.
//' @param gcv_identification_check boolean `TRUE` for a loose check for model
//' identification; `FALSE` otherwise. Default is `FALSE`.
//' @param target_level_kp2006_rank_test (only relevant if `gcv_identification_check` is
//' `TRUE`) numeric level of the Kleibergen Paap 2006 rank test. If it is
//' strictly grater than zero, then the iterative Kleibergen Paap 2006 rank
//' test at `level = target_level_kp2006_rank_test / n_factors` is used to compute an initial estimator
//' of the rank of the factor loadings in the Chen Fang 2019 rank test.
//' Otherwise, the initial rank estimator is taken to be the number of singular
//' values above `n_observations^(-1/4)`. Default is `0.05` (as correction
//' for multiple testing).
//' @param relaxed boolean `TRUE` if you want to compute a post-selection
//' unpenalized tradable factor risk premia to remove the bias due to shrinkage;
//' FALSE` otherwise. Default is `FALSE`.
//' @param include_standard_errors boolean `TRUE` if you want to compute the
//' oracle factor risk premia HAC standard errors; `FALSE` otherwise.
//' Default is `FALSE`.
//'
//' @return a list containing the `n_factors`-dimensional vector of oracle
//' tradable factor risk premia in `"risk_premia"`, and the optimal penalty
//' parameter value in `"penalty_parameter"`.
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OracleTFRPGCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const bool one_stddev_rule = false,
  const bool gcv_scaling_n_assets = false,
  const bool gcv_identification_check = false,
  const double target_level_kp2006_rank_test = 0.05,
  const bool relaxed = false,
  const bool include_standard_errors = false
);

//' Compute optimal oracle tradable factor risk premia under cross validation
//'
//' @name OracleTFRPCVCpp
//' @description Computes optimal oracle tradable factor risk premia based
//' on moments extracted from factors and test asset excess returns and oracle
//' weights over various penalty parameter values. Tuning is performed via
//' Cross Validation (CV). Oracle weights can be based on the correlation
//' between factors and returns, on the regression coefficients of returns on
//' factors or on the first-step tradable risk premia estimator.
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
//' @param weighting_type character specifying the type of oracle weights:
//' based on the correlation between factors and returns `'c'`; based on the
//' regression coefficients of returns on factors `'b'`; based on the first-step
//' tradable risk premia estimator `'a'`; otherwise a vector of ones (any other
//' character). Default is `'c'`.
//' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
//' whose score is not higher than one standard error above the score of the
//' best model; `FALSE` for picking the best model. Default is `FALSE`.
//' @param n_folds integer number of k-fold for cross validation. Default is `5`.
//' @param relaxed boolean `TRUE` if you want to compute a post-selection
//' unpenalized tradable factor risk premia to remove the bias due to shrinkage;
//' FALSE` otherwise. Default is `FALSE`.
//' @param include_standard_errors boolean `TRUE` if you want to compute the
//' oracle factor risk premia HAC standard errors; `FALSE` otherwise.
//' Default is `FALSE`.
//'
//' @return a list containing the n_factors-dimensional vector of oracle
//' tradable factor risk premia in "risk_premia", and the optimal penalty
//' parameter value in "penalty_parameter".
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OracleTFRPCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const bool one_stddev_rule = false,
  const unsigned int n_folds = 5,
  const bool relaxed = false,
  const bool include_standard_errors = false
);

//' Compute optimal oracle tradable factor risk premia under rolling validation
//'
//' @name OracleTFRPRVCpp
//' @description Computes optimal oracle tradable factor risk premia based
//' on moments extracted from factors and test asset excess returns and oracle
//' weights over various penalty parameter values. Tuning is performed via
//' Rolling Validation (RV). Oracle weights can be based on the correlation
//' between factors and returns, on the regression coefficients of returns on
//' factors or on the first-step tradable risk premia estimator.
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
//' @param weighting_type character specifying the type of oracle weights:
//' based on the correlation between factors and returns `'c'`; based on the
//' regression coefficients of returns on factors `'b'`; based on the first-step
//' tradable risk premia estimator `'a'`; otherwise a vector of ones (any other
//' character). Default is `'c'`.
//' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
//' whose score is not higher than one standard error above the score of the
//' best model; `FALSE` for picking the best model. Default is `FALSE`.
//' @param n_train_observations number of observations in the rolling training
//' set. Default is `120`.
//' @param n_test_observations number of observations in the test set. Default
//' is `12`.
//' @param roll_shift number of observation shift when moving from the rolling
//' window to the next one. Default is `12`.
//' @param relaxed boolean `TRUE` if you want to compute a post-selection
//' unpenalized tradable factor risk premia to remove the bias due to shrinkage;
//' FALSE` otherwise. Default is `FALSE`.
//' @param include_standard_errors boolean `TRUE` if you want to compute the
//' oracle factor risk premia HAC standard errors; `FALSE` otherwise.
//' Default is `FALSE`.
//'
//' @return a list containing the n_factors-dimensional vector of oracle
//' tradable factor risk premia in "risk_premia", and the optimal penalty
//' parameter value in "penalty_parameter".
//'
//' @noRd
//'
// [[Rcpp::export]]
Rcpp::List OracleTFRPRVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const bool one_stddev_rule = false,
  const unsigned int n_train_observations = 120,
  const unsigned int n_test_observations = 12,
  const unsigned int roll_shift = 12,
  const bool relaxed = false,
  const bool include_standard_errors = false
);

//' Compute oracle tradable factor risk premia
//'
//' @name OracleTFRPCpp
//' @description Computes oracle tradable factor risk premia based on
//' pre-computed tradable factor risk premia and oracle penalty weights for
//' various penalty parameter values.
//'
//' @param tradable_frp `n_factors`-dimensional vector of tradable factor risk premia.
//' Usually, it is the output of function `TFRPCpp`
//' @param weights `n_factors`-dimensional vector of weights for the penalty
//' term. Usually, it is the output of function `OracleWeightsCpp`.
//' @param penalty_parameters `n_parameters`-dimensional vector of penalty
//' parameter values from smallest to largest.
//'
//' @return `n_factors x n_parameters`-dimensional matrix of oracle
//' tradable factor risk premia.
//'
//' @noRd
arma::mat OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const arma::vec& penalty_parameter
);

// For internal use
// Computes oracle tradable factor risk premia based on pre-computed
// tradable factor risk premia and oracle penalty weights for a specific
// penalty parameter value.
arma::vec OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const double penalty_parameter
);

// Function for internal use
// Computes the Relaxed Oracle tradable factor risk premia.
arma::vec RelaxOracleTFRP(
  const arma::uvec& idx_nonzero,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

//' Compute the HAC standard errors of the nonzero oracle tradable factor risk
//' premia
//'
//' @name StandardErrorsOracleTFRPCpp
//' @description Computes the HAC standard errors of oracle tradable factor
//' risk premia based on moments extracted from factors and test asset excess
//' returns. It uses the Newey-West (1994) plug-in procedure to select the
//' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//'
//' @param idx_nonzero index vector of nonzero oracle tradable factor risk premia.
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
//'
//' @noRd
//'
//' @return `n_factors`-dimensional vector of standard errors for oracle
//' tradable factor risk premia.
arma::vec StandardErrorsOracleTFRPCpp(
  const arma::uvec& idx_nonzero,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// Function for internal use
// Computes the optimal tuning parameter according to the one "standard deviation
// rule".
unsigned int ComputeOneStdDevRuleCpp(
    const arma::vec& score
);

#endif
