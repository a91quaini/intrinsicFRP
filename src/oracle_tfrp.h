// Author: Alberto Quaini

#ifndef ORACLE_TFRP_H
#define ORACLE_TFRP_H

#include <RcppArmadillo.h>

//////////////////////////////////////////////////////////
////////  Oracle Tradable Factor Risk Premia /////////////
//////////////////////////////////////////////////////////

// @title Oracle tradable factor risk premia.
//
// @name OracleTFRPCpp
// @description Computes Oracle tradable factor risk premia of
// Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683> from data on
// `K` factors `F = [F_1,...,F_K]'` and test asset excess returns `R`:
// `OTFRP = argmin_x ||TFRP - x||_2^2 + tau * sum_{k=1}^K w_k * |x_k|`,
// where `TFRP` is the tradable factor risk premia estimator, `tau > 0` is a
// penalty parameter, and the Oracle weights are given by
// `w_k = 1 / ||corr[F_k, R]||_2^2`.
// This estimator is called "Oracle" in the sense that the probability that
// the index set of its nonzero estimated risk premia equals the index set of
// the true strong factors tends to 1 (Oracle selection), and that on the strong
// factors, the estimator reaches the optimal asymptotic Normal distribution.
// Here, strong factors are those that have a nonzero population marginal
// correlation with asset excess returns.
// Tuning of the penalty parameter `tau` is performed via Generalized Cross
// Validation (GCV), Cross Validation (CV) or Rolling Validation (RV).
// GCV tunes parameter `tau` by minimizing the criterium:
// `||PE(tau)||_2^2 / (1-df(tau)/T)^2`
// where
// `PE(tau) = E[R] - beta_{S(tau)} * OTFRP(tau)`
// are the pricing errors of the model for given tuning parameter `tau`,
// with `S(tau)` being the index set of the nonzero Oracle TFRP computed with
// tuning parameter `tau`, and
// `beta_{S(tau)} = Cov[R, F_{S(tau)}] * (Cov[F_{S(tau)}, R] * V[R]^{-1} * Cov[R, F_{S(tau)}])^{-1}`
// the regression coefficients of the test assets excess returns on the
// factor mimicking portfolios,
// and `df(tau) = |S(tau)|` are the degrees of freedom of the model, given by the
// number of nonzero Oracle TFRP.
// CV and RV, instead, choose the value of `tau` that minimize the criterium:
// `PE(tau)' * V[PE(tau)]^{-1} PE(tau)`
// where `V[PE(tau)]` is the diagonal matrix collecting the marginal variances
// of pricing errors `PE(tau)`, and each of these components are
// aggregated over k-fold cross-validated data or over rolling windows of data,
// respectively.
// Oracle weights can be based on the correlation between factors and returns
// (suggested approach),
// on the regression coefficients of returns on factors or on the first-step
// tradable risk premia estimator. Optionally computes the corresponding
// heteroskedasticity and autocorrelation robust standard errors using the
// Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the number
// of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
// For the standard error computations, the function allows to internally
// pre-whiten the series by fitting a VAR(1),
// i.e., a vector autoregressive model of order 1.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>.
//
// @param returns A `n_observations x n_returns`-dimensional matrix of test asset
// excess returns.
// @param factors A `n_observations x n_factors`-dimensional matrix of factors.
// @param penalty_parameters A `n_parameters`-dimensional vector of penalty
// parameter values from smallest to largest.
// @param weighting_type A character specifying the type of adaptive weights:
// based on the correlation between factors and returns `'c'`; based on the
// regression coefficients of returns on factors `'b'`; based on the first-step
// tradable risk premia estimator `'a'`; otherwise a vector of ones (any other
// character). Default is `'c'`.
// @param tuning_type A character indicating the parameter tuning type:
// `'g'` for generalized cross validation; `'c'` for K-fold cross validation;
// `'r'` for rolling validation. Default is `'g'`.
// @param one_stddev_rule A boolean: `TRUE` for picking the most parsimonious model
// whose score is not higher than one standard error above the score of the
// best model; `FALSE` for picking the best model. Default is `TRUE`.
// @param gcv_scaling_n_assets (only relevant for `tuning_type ='g'`)
// A boolean: `TRUE` for sqrt(n_assets) scaling (`sqrt(n_assets) / n_observations`);
// `FALSE` otherwise (`1 / n_observations`). Default is `FALSE`.
// @param gcv_identification_check (only relevant for `tuning_type ='g'`)
// A boolean: `TRUE` for a loose check for model
// identification; `FALSE` otherwise. Default is `FALSE`.
// @param target_level_kp2006_rank_test (only relevant for `tuning_type ='g'`
// and if `gcv_identification_check` is
// `TRUE`) A number indicating the level of the Kleibergen Paap 2006 rank test.
// If it is strictly grater than zero, then the iterative Kleibergen Paap 2006 rank
// test at `level = target_level_kp2006_rank_test / n_factors`
// (where division by the number of factors is performed as a Bonferroni correction for
// multiple testing) is used to compute an initial estimator
// of the rank of the factor loadings in the Chen Fang 2019 rank test.
// Otherwise, the initial rank estimator is taken to be the number of singular
// values above `n_observations^(-1/4)`.
// Default is `0.05`.
// @param n_folds (only relevant for `tuning_type ='c'`) An integer indicating
// the number of k-fold for cross validation. Default is `5`.
// @param n_train_observations (only relevant for `tuning_type ='r'`) The number of
// observations in the rolling training set. Default is `120`.
// @param n_test_observations (only relevant for `tuning_type ='r'`) The number of
// observations in the test set. Default is `12`.
// @param roll_shift (only relevant for `tuning_type ='r'`) The number of observation
// shift when moving from the rolling window to the next one. Default is `12`.
// @param relaxed A boolean: `TRUE` if you want to compute a post-selection
// unpenalized tradable factor risk premia to remove the bias due to shrinkage;
// FALSE` otherwise. Default is `FALSE`.
// @param include_standard_errors A boolean `TRUE` if you want to compute the
// adaptive tradable factor risk premia HAC standard errors; `FALSE` otherwise.
// Default is `FALSE`.
// @param hac_prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// robust covariance (HAC) estimation. Default is `false`.
// @param plot_score A boolean: `TRUE` for plotting the model score; `FALSE` otherwise.
// Default is `TRUE`.
// @param check_arguments A boolean `TRUE` for internal check of all function
// arguments; `FALSE` otherwise. Default is `TRUE`.
//
// @return A list containing the `n_factors`-dimensional vector of adaptive
// tradable factor risk premia in `"risk_premia"`; the optimal penalty
// parameter value in `"penalty_parameter"`; the model score for each penalty
// parameter value in `"model_score"`;  if `include_standard_errors = TRUE`, then
// it further includes `n_factors`-dimensional vector of tradable factor risk
// premia standard errors in `"standard_errors"`.
//
// [[Rcpp::export]]
Rcpp::List OracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const char tuning_type,
  const bool one_stddev_rule,
  const bool gcv_scaling_n_assets,
  const bool gcv_identification_check,
  const double target_level_kp2006_rank_test,
  const unsigned int n_folds,
  const unsigned int n_train_observations,
  const unsigned int n_test_observations,
  const unsigned int roll_shift,
  const bool relaxed,
  const bool include_standard_errors,
  const bool hac_prewhite
);

// Function for internal use.
// Computes optimal adaptive tradable factor risk premia
// using generalized cross validation parameter tuning.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
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
  const bool include_standard_errors = false,
  const bool hac_prewhite = false
);

// Function for internal use.
// Computes optimal adaptive tradable factor risk premia
// using cross validation parameter tuning.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
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
  const bool include_standard_errors = false,
  const bool hac_prewhite = false
);

// Function for internal use.
// Computes optimal adaptive tradable factor risk premia
// using rolling validation parameter tuning.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
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
  const bool include_standard_errors = false,
  const bool hac_prewhite = false
);

// Function for internal use
// Computes oracle tradable factor risk premia
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
arma::mat OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const arma::vec& penalty_parameter
);

// For internal use
// Computes oracle tradable factor risk premia based on pre-computed
// tradable factor risk premia and oracle penalty weights for a specific
// penalty parameter value.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
arma::vec OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const double penalty_parameter
);

// Function for internal use
// Computes the Relaxed Oracle tradable factor risk premia.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
arma::vec RelaxOracleTFRP(
  const arma::uvec& idx_nonzero,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// Compute the HAC standard errors of the nonzero oracle tradable factor risk
// premia
//
// @name StandardErrorsOracleTFRPCpp
// @description Computes the HAC standard errors of oracle tradable factor
// risk premia based on moments extracted from factors and test asset excess
// returns. It uses the Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the
// number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
//
// @param idx_nonzero index vector of nonzero oracle tradable factor risk premia.
// @param returns `n_observations x n_returns`-dimensional matrix of test asset
// excess returns.
// @param factors `n_observations x n_factors`-dimensional matrix of risk
// factors.
// @param covariance_factors_returns `n_factors x n_returns`-dimensional
// covariance matrix between factors and test asset excess returns.
// @param variance_returns `n_returns x n_returns`-dimensional covariance
// matrix of test asset excess returns.
// @param mean_returns `n_returns`-dimensional mean vector of test asset excess
// returns.
// @param hac_prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// robust covariance (HAC) estimation. Default is `false`.
//
// @return `n_factors`-dimensional vector of standard errors for oracle
// tradable factor risk premia.
//
arma::vec StandardErrorsOracleTFRPCpp(
  const arma::uvec& idx_nonzero,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite = false
);

// Function for internal use
// Computes the optimal tuning parameter according to the one "standard deviation
// rule".
// All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>
unsigned int ComputeOneStdDevRuleCpp(
  const arma::vec& score
);

#endif
