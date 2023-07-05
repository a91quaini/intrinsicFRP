# Author: Alberto Quaini

##################################################
####  Adaptive Intrinsic Factor Risk Premia ######
##################################################

#' Compute optimal adaptive intrinsic factor risk premia
#'
#' @name OptimalAdaptiveIFRP
#' @description Computes optimal adaptive intrinsic factor risk premia based on
#' pre-computed intrinsic factor risk premia for
#' various penalty parameter values. Tuning is performed via
#' Generalized Cross Validation (GCV), Cross Validation (CV) or Rolling
#' Validation (RV). Adaptive weights can be based on the correlation between
#' factors and returns, on the regression coefficients of returns on factors
#' or on the first-step intrinsic risk premia estimator. Optionally computes
#' the corresponding heteroskedasticity and autocorrelation robust standard
#' errors using the Newey-West (1994) plug-in procedure to select the number
#' of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of factors.
#' @param penalty_parameters `n_parameters`-dimensional vector of penalty
#' parameter values from smallest to largest.
#' @param weighting_type character specifying the type of adaptive weights:
#' based on the correlation between factors and returns `'c'`; based on the
#' regression coefficients of returns on factors `'b'`; based on the first-step
#' intrinsic risk premia estimator `'a'`; otherwise a vector of ones (any other
#' character). Default is `'c'`.
#' @param tuning_type character indicating the parameter tuning type: `'g'` for
#' generalized cross validation; `'r'` for rolling validation. Default is `'g'`.
#' @param include_standard_errors boolean `TRUE` if you want to compute the
#' adaptive intrinsic factor risk premia HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param one_stddev_rule boolean `TRUE` for picking the most parsimonious model
#' whose score is not higher than one standard error above the score of the
#' best model; `FALSE` for picking the best model. Default is `FALSE`.
#' @param gcv_vr_weighting boolean `TRUE` for scaling pricing errors by
#' the inverse variance matrix of asset excess returns; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param gcv_scaling_n_assets (only relevant for `tuning_type ='g'`)
#' boolean `TRUE` for sqrt(n_assets) scaling (`sqrt(n_assets) / n_observations`);
#' `FALSE` otherwise (`1 / n_observations`). Default is `FALSE`.
#' @param gcv_identification_check (only relevant for `tuning_type ='g'`)
#' boolean `TRUE` for a loose check for model identification; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param target_level_kp2006_rank_test (only relevant for `tuning_type ='g'` if
#' `gcv_identification_check` is `TRUE`) numeric level of the Kleibergen Paap
#' 2006 rank test. If it is strictly grater than zero, then the iterative
#' Kleibergen Paap 2006 rank test at
#' `level = target_level_kp2006_rank_test / n_factors` (as a correction for multiple testing) is used to
#' compute an initial estimator of the rank of the factor loadings in the Chen
#' Fang 2019 rank test. Otherwise, the initial rank estimator is taken to be
#' the number of singular values above `n_observations^(-1/4)`. Default is
#' `0.05` (as correction for multiple testing).
#' `gcv_identification_check` is `TRUE`) numeric level of the Kleibergen
#' Paap 2006 rank test. If it is strictly grater than zero, then the iterative
#' Kleibergen Paap 2006 rank test at level `level_kp2005_test` is used to
#' compute an initial estimator of the rank of the factor loadings in the
#' Chen Fang 2019 rank test. Otherwise, the initial rank estimator is taken
#' to be the number of singular values above `n_observations^(-1/3)`. Default
#' is `0.005` (as correction for multiple testing).
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
#' intrinsic factor risk premia in `"risk_premia"`; the optimal penalty
#' parameter value in `"penalty_parameter"`; the model score for each penalty
#' parameter value in `"model_score"`;  if `include_standard_errors=TRUE`, then
#' it further includes `n_factors`-dimensional vector of intrinsic factor risk
#' premia standard errors in `"standard_errors"`.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' penalty_parameters = seq(0., 1., length.out = 100)
#'
#' # compute optimal adaptive intrinsic factor risk premia and their standard errors
#' aifrp = OptimalAdaptiveIFRP(
#' returns,
#' factors,
#' penalty_parameters,
#' include_standard_errors = TRUE
#' )
#'
#' @export
OptimalAdaptiveIFRP = function(
  returns,
  factors,
  penalty_parameters,
  weighting_type = 'c',
  tuning_type = 'g',
  include_standard_errors = FALSE,
  one_stddev_rule = FALSE,
  gcv_vr_weighting = FALSE,
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

  stopifnot("`check_arguments` must be boolean" = is.logical(check_arguments))
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`penalty_parameters` contains non-numeric values" = is.numeric(penalty_parameters))
    stopifnot("`penalty_parameters` contains missing values (NA/NaN)" = !anyNA(penalty_parameters))
    stopifnot("`weighting_type` must be a character value" = is.character(weighting_type))
    stopifnot("`tuning_type` must be a character value" = is.character(tuning_type))
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))
    stopifnot("`one_stddev_rule` must be boolean" = is.logical(one_stddev_rule))
    stopifnot("`gcv_vr_weighting` must be boolean" = is.logical(gcv_vr_weighting))
    stopifnot("`gcv_scaling_n_assets` must be boolean" = is.logical(gcv_scaling_n_assets))
    stopifnot("`gcv_identification_check` must be boolean" = is.logical(gcv_identification_check))
    stopifnot("`target_level_kp2006_rank_test` contains non-numeric values" = is.numeric(target_level_kp2006_rank_test))
    stopifnot("`n_folds` contains non-numeric values" = is.numeric(n_folds))
    stopifnot("`n_folds` should be between 2 and n_returns" = n_folds > 2 || n_folds < nrow(returns))
    stopifnot("`n_train_observations` contains non-numeric values" = is.numeric(n_train_observations))
    stopifnot("`n_train_observations` should be between 10 and n_obervations - n_test_observations" = n_train_observations > 10 || n_train_observations < nrow(returns) - n_test_observations)
    stopifnot("`n_test_observations` contains non-numeric values" = is.numeric(n_test_observations))
    stopifnot("`n_test_observations` should be between 10 and n_observations/2" = n_test_observations > 10 || n_test_observations < nrow(returns) / 2)
    stopifnot("`roll_shift` contains non-numeric values" = is.numeric(roll_shift))
    stopifnot("`roll_shift` should be between 1 and n_test_observations" = roll_shift >= 1 || roll_shift < n_test_observations)
    stopifnot("`plot_score` must be boolean" = is.logical(plot_score))
    penalty_parameters = sort(penalty_parameters)

  }

  output = switch(
    tuning_type,
    'g' = {

      .Call(`_intrinsicFRP_OptimalAdaptiveIFRPGCVCpp`,
        returns,
        factors,
        stats::cov(factors, returns),
        stats::cov(returns),
        colMeans(returns),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        gcv_vr_weighting,
        gcv_scaling_n_assets,
        gcv_identification_check,
        target_level_kp2006_rank_test,
        include_standard_errors
      )

    },
    'c' = {

      .Call(`_intrinsicFRP_OptimalAdaptiveIFRPCVCpp`,
        returns,
        factors,
        stats::cov(factors, returns),
        stats::cov(returns),
        colMeans(returns),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        n_folds,
        include_standard_errors
      )

    },
    'r' = {

      .Call(`_intrinsicFRP_OptimalAdaptiveIFRPRVCpp`,
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
        include_standard_errors
      )

    },
    stop("Invalid `tuning_type` value")
  )

  if (plot_score) {PlotAdaptiveIFRPModelScore(output, penalty_parameters)}

  return(output)

}

#' Compute adaptive intrinsic factor risk premia
#'
#' @name AdaptiveIFRP
#' @description Computes adaptive intrinsic factor risk premia with user-defined
#' weights for various penalty parameter values.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of factors.
#' @param penalty_parameters `n_parameters`-dimensional vector of penalty
#' parameter values from smallest to largest.
#' @param weights `n_factors`-dimensional vector of weights determining a
#' separate penalty parameter for each risk premium, given by vector
#' `penalty_parameter * weights`.
#' Default is a vector of ones, i.e., the same penalty parameter
#' `penalty_parameter` is applied to each risk premium.
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return `n_factors x n_parameters`-dimensional matrix of adaptive
#' intrinsic factor risk premia.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' weights = 1. / rowSums(stats::cor(factors, returns))
#' penalty_parameters = seq(0., 1., length.out = 100)
#'
#' # compute adaptive intrinsic factor risk premia for each penalty parameter
#' aifrp = AdaptiveIFRP(
#' returns,
#' factors,
#' penalty_parameters,
#' weights
#' )
#'
#' @export
AdaptiveIFRP = function(
  returns,
  factors,
  penalty_parameters,
  weights = rep(1., ncol(factors)),
  check_arguments = TRUE
) {

  stopifnot("`check_arguments` must be boolean" = is.logical(check_arguments))
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`penalty_parameters` contains non-numeric values" = is.numeric(penalty_parameters))
    stopifnot("`penalty_parameters` contains missing values (NA/NaN)" = !anyNA(penalty_parameters))
    stopifnot("`weights` contains non-numeric values" = is.numeric(weights))
    stopifnot("`weights` contains missing values (NA/NaN)" = !anyNA(weights))

  }

  ifrp = .Call(`_intrinsicFRP_IFRPCpp`,
    returns,
    factors,
    FALSE
  )

  return(.Call(`_intrinsicFRP_AdaptiveIFRPCpp`,
    ifrp$risk_premia,
    weights,
    penalty_parameters
  ))

}
