# Author: Alberto Quaini

#######################################
######  Intrinsic Factor Risk Premia ##
#######################################

#' Compute intrinsic factor risk premia
#'
#' @name IFRP
#' @description Computes intrinsic factor risk premia from data on factors and
#' test asset excess returns. Optionally computes the corresponding
#' heteroskedasticity and autocorrelation robust standard errors using the
#' Newey-West (1994) plug-in procedure to select the number of relevant lags,
#' i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of factors.
#' @param include_standard_errors boolean `TRUE` if you want to compute the
#' intrinsic factor risk premia HAC standard errors; `FALSE` otherwise. Default
#' is `FALSE`.
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing `n_factors`-dimensional vector of intrinsic factor
#' risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`, then
#' it further includes `n_factors`-dimensional vector of intrinsic factor risk
#' premia standard errors in `"standard_errors"`.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute intrinsic factor risk premia and their standard errors
#' ifrp = IFRP(returns, factors, include_standard_errors = TRUE)
#'
#' @export
IFRP = function(
  returns,
  factors,
  include_standard_errors = FALSE,
  check_arguments = TRUE
) {

  stopifnot("`check_arguments` must be boolean" = is.logical(check_arguments))
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))

  }

  if (include_standard_errors) {

    covariance_factors_returns = stats::cov(factors, returns)
    variance_returns = stats::cov(returns)
    mean_returns = colMeans(returns)

    ifrp = .Call(`_intrinsicFRP_IFRPCpp`,
      covariance_factors_returns,
      variance_returns,
      mean_returns
    )

    return(list(
      "risk_premia" = ifrp,
      "standard_errors" = .Call(`_intrinsicFRP_StandardErrorsIFRPCpp`,
        ifrp,
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        colMeans(factors)
      )
    ))

  }

  return(list(
    "risk_premia" = .Call(`_intrinsicFRP_IFRPCpp`,
      stats::cov(factors, returns),
      stats::cov(returns),
      colMeans(returns)
  )))

}

