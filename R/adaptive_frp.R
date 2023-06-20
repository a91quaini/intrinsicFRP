# Author: Alberto Quaini

########################################
####  Adaptive Factor Risk Premia ######
########################################

#' Compute optimal adaptive factor risk premia
#'
#' @name OptimalAdaptiveFRP
#' @description Computes optimal adaptive factor risk premia based on
#' pre-computed intrinsic factor risk premia for
#' various penalty parameter values. Tuning is performed via
#' minimizing the ratio between model misspecification and mis-identification.
#' Adaptive weights can be based on the correlation between
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
#' @param include_standard_errors boolean `TRUE` if you want to compute the
#' adaptive intrinsic factor risk premia HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param plot_score boolean `TRUE` for plotting the model score; `FALSE` otherwise.
#' Default is `TRUE`.
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing the `n_factors`-dimensional vector of adaptive
#'factor risk premia in `"risk_premia"`; the optimal penalty
#' parameter value in `"penalty_parameter"`; the model score for each penalty
#' parameter value in `"model_score"`;  if `include_standard_errors=TRUE`, then
#' it further includes `n_factors`-dimensional vector of factor risk
#' premia standard errors in `"standard_errors"`.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' penalty_parameters = seq(0., 1., length.out = 100)
#'
#' # compute optimal adaptive factor risk premia and their standard errors
#' aifrp = OptimalAdaptiveFRP(
#' returns,
#' factors,
#' penalty_parameters,
#' include_standard_errors = FALSE
#' )
#'
#' @export
OptimalAdaptiveFRP = function(
  returns,
  factors,
  penalty_parameters,
  weighting_type = 'c',
  include_standard_errors = FALSE,
  plot_score = TRUE,
  check_arguments = TRUE
) {

  stopifnot("`check_arguments` must be boolean" = is.logical(check_arguments))
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`penalty_parameters` contains non-numeric values" = is.numeric(penalty_parameters))
    stopifnot("`penalty_parameters` contains missing values (NA/NaN)" = !anyNA(penalty_parameters))
    stopifnot("`weighting_type` must be a character value" = is.character(weighting_type))
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))
    stopifnot("`plot_score` must be boolean" = is.logical(plot_score))
    penalty_parameters = sort(penalty_parameters)

  }

  covariance_factors_returns = stats::cov(factors, returns)
  variance_returns = stats::cov(returns)
  mean_returns = colMeans(returns)

  output = .Call(`_intrinsicFRP_OptimalAdaptiveFRPCpp`,
    returns,
    factors,
    covariance_factors_returns,
    variance_returns,
    mean_returns,
    penalty_parameters,
    weighting_type
  )

  if (include_standard_errors) {

    output[["standard_errors"]] = .Call(`_intrinsicFRP_StandardErrorsAdaptiveFRPCpp`,
      output$risk_premia,
      returns,
      factors,
      variance_returns,
      mean_returns
    )

  }

  if (plot_score) {PlotAdaptiveIFRPModelScore(output, penalty_parameters)}

  return(output)

}
