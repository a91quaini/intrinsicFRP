# Author: Alberto Quaini

##################################################
####  Oracle Tradable Factor Risk Premia ######
##################################################

#' Compute optimal adaptive tradable factor risk premia
#'
#' @name OracleTFRP
#' @description Computes optimal adaptive tradable factor risk premia for
#' various penalty parameter values from data on factors and test asset excess
#' returns. Tuning is performed via Generalized Cross Validation (GCV),
#' Cross Validation (CV) or Rolling Validation (RV). Oracle weights can be based
#' on the correlation between factors and returns, on the regression
#' coefficients of returns on factors or on the first-step tradable risk premia
#' estimator. Optionally computes the corresponding heteroskedasticity and
#' autocorrelation robust standard errors using the Newey-West (1994) plug-in
#' procedure to select the number of relevant lags, i.e.,
#' `n_lags = 4 * (n_observations/100)^(2/9)`.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of factors.
#' @param penalty_parameters `n_parameters`-dimensional vector of penalty
#' parameter values from smallest to largest.
#' @param weighting_type character specifying the type of adaptive weights:
#' based on the correlation between factors and returns `'c'`; based on the
#' regression coefficients of returns on factors `'b'`; based on the first-step
#' tradable risk premia estimator `'a'`; otherwise a vector of ones (any other
#' character). Default is `'c'`.
#' @param tuning_type character indicating the parameter tuning type: `'g'` for
#' generalized cross validation; `'r'` for rolling validation. Default is `'g'`.
#' @param include_standard_errors boolean `TRUE` if you want to compute the
#' adaptive tradable factor risk premia HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param relaxed boolean `TRUE` if you want to compute a post-selection
#' unpenalized tradable factor risk premia to remove the bias due to shrinkage;
#' FALSE` otherwise. Default is `FALSE`.
#' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
#' whose score is not higher than one standard error above the score of the
#' best model; `FALSE` for picking the best model. Default is `FALSE`.
#' @param gcv_scaling_n_assets (only relevant for `tuning_type ='g'`)
#' boolean `TRUE` for sqrt(n_assets) scaling (`sqrt(n_assets) / n_observations`);
#' `FALSE` otherwise (`1 / n_observations`). Default is `FALSE`.
#' @param gcv_identification_check (only relevant for `tuning_type ='g'`)
#' boolean `TRUE` for a loose check for model
#' identification; `FALSE` otherwise. Default is `FALSE`.
#' @param target_level_kp2006_rank_test (only relevant for `tuning_type ='g'`
#' and if `gcv_identification_check` is
#' `TRUE`) numeric level of the Kleibergen Paap 2006 rank test. If it is
#' strictly grater than zero, then the iterative Kleibergen Paap 2006 rank
#' test at `level = target_level_kp2006_rank_test / n_factors` is used to compute an initial estimator
#' of the rank of the factor loadings in the Chen Fang 2019 rank test.
#' Otherwise, the initial rank estimator is taken to be the number of singular
#' values above `n_observations^(-1/4)`. Default is `0.05` (as correction
#' for multiple testing).
#' @param n_folds (only relevant for `tuning_type ='c'`) integer number of k-fold
#' for cross validation. Default is `5`.
#' @param n_train_observations (only relevant for `tuning_type ='r'`) number of
#' observations in the rolling training set. Default is `120`.
#' @param n_test_observations (only relevant for `tuning_type ='r'`) number of
#' observations in the test set. Default is `12`.
#' @param roll_shift (only relevant for `tuning_type ='r'`) number of observation
#' shift when moving from the rolling window to the next one. Default is `12`.
#' @param plot_score boolean `TRUE` for plotting the model score; `FALSE` otherwise.
#' Default is `TRUE`.
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing the `n_factors`-dimensional vector of adaptive
#' tradable factor risk premia in `"risk_premia"`; the optimal penalty
#' parameter value in `"penalty_parameter"`; the model score for each penalty
#' parameter value in `"model_score"`;  if `include_standard_errors=TRUE`, then
#' it further includes `n_factors`-dimensional vector of tradable factor risk
#' premia standard errors in `"standard_errors"`.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' penalty_parameters = seq(0., 1., length.out = 100)
#'
#' # compute optimal adaptive tradable factor risk premia and their standard errors
#' oracle_tfrp = OracleTFRP(
#' returns,
#' factors,
#' penalty_parameters,
#' include_standard_errors = TRUE
#' )
#'
#' @export
OracleTFRP = function(
  returns,
  factors,
  penalty_parameters,
  weighting_type = 'c',
  tuning_type = 'g',
  relaxed = FALSE,
  include_standard_errors = FALSE,
  one_stddev_rule = FALSE,
  gcv_scaling_n_assets = FALSE,
  gcv_identification_check = FALSE,
  target_level_kp2006_rank_test = 0.05,
  n_folds = 5,
  n_train_observations = 120,
  n_test_observations = 12,
  roll_shift = 12,
  plot_score = TRUE,
  check_arguments = TRUE
) {

  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`penalty_parameters` contains non-numeric values" = is.numeric(penalty_parameters))
    stopifnot("`penalty_parameters` contains missing values (NA/NaN)" = !anyNA(penalty_parameters))
    stopifnot("`weighting_type` must be a character value" = is.character(weighting_type))
    stopifnot("`tuning_type` must be a character value" = is.character(tuning_type))
    stopifnot("`relaxed` must be boolean" = is.logical(relaxed))
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))
    stopifnot("`one_stddev_rule` must be boolean" = is.logical(one_stddev_rule))
    stopifnot("`gcv_scaling_n_assets` must be boolean" = is.logical(gcv_scaling_n_assets))
    stopifnot("`gcv_identification_check` must be boolean" = is.logical(gcv_identification_check))
    stopifnot("`target_level_kp2006_rank_test` must be numeric" = is.numeric(target_level_kp2006_rank_test))
    stopifnot("`target_level_kp2006_rank_test` must be between 0 and 1" = (target_level_kp2006_rank_test >= 0.) & (target_level_kp2006_rank_test <= 1.))
    stopifnot("`n_folds` must be numeric" = is.numeric(n_folds))
    stopifnot("`n_folds` should be between 2 and n_returns" = n_folds > 2 || n_folds < nrow(returns))
    stopifnot("`n_train_observations` must be numeric" = is.numeric(n_train_observations))
    stopifnot("`n_train_observations` should be between 10 and n_obervations - n_test_observations" = n_train_observations > 10 || n_train_observations < nrow(returns) - n_test_observations)
    stopifnot("`n_test_observations` must be numeric" = is.numeric(n_test_observations))
    stopifnot("`n_test_observations` should be between 10 and n_observations/2" = n_test_observations > 10 || n_test_observations < nrow(returns) / 2)
    stopifnot("`roll_shift` must be numeric" = is.numeric(roll_shift))
    stopifnot("`roll_shift` should be between 1 and n_test_observations" = roll_shift >= 1 || roll_shift < n_test_observations)
    stopifnot("`plot_score` must be boolean" = is.logical(plot_score))
    penalty_parameters = sort(penalty_parameters)

  }

  output = switch(
    tuning_type,
    'g' = {

      .Call(`_intrinsicFRP_OracleTFRPGCVCpp`,
        returns,
        factors,
        stats::cov(factors, returns),
        stats::cov(returns),
        colMeans(returns),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        gcv_scaling_n_assets,
        gcv_identification_check,
        target_level_kp2006_rank_test,
        relaxed,
        include_standard_errors
      )

    },
    'c' = {

      .Call(`_intrinsicFRP_OracleTFRPCVCpp`,
        returns,
        factors,
        stats::cov(factors, returns),
        stats::cov(returns),
        colMeans(returns),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        n_folds,
        relaxed,
        include_standard_errors
      )

    },
    'r' = {

      .Call(`_intrinsicFRP_OracleTFRPRVCpp`,
        returns,
        factors,
        stats::cov(factors, returns),
        stats::cov(returns),
        colMeans(returns),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        n_train_observations,
        n_test_observations,
        roll_shift,
        relaxed,
        include_standard_errors
      )

    },
    stop("Invalid `tuning_type` value")
  )

  if (plot_score) {PlotOracleTFRPModelScore(output, penalty_parameters)}

  return(output)

}
