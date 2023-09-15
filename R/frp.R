# Author: Alberto Quaini

###################################
######  Factor Risk Premia ########
###################################

#' Compute factor risk premia from data
#'
#' @name FRP
#' @description Computes Fama MachBeth (1973) or misspecification-robust factor
#' risk premia of Kan Robotti Shanken (2013) from data on factors and test
#' asset excess returns. Optionally computes the corresponding
#' heteroskedasticity and autocorrelation robust standard errors using the
#' Newey-West (1994) plug-in procedure to select the number of relevant lags,
#' i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of factors.
#' @param misspecification_robust boolean `TRUE` for the
#' "misspecification-robust" Kar Robotti Shanken GLS approach using the
#' inverse covariance matrix of returns; `FALSE` for standard Fama-Mac-Beth
#' risk premia. Default is `TRUE`.
#' @param include_standard_errors boolean `TRUE` if you want to compute the
#' factor risk premia HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing `n_factors`-dimensional vector of factor
#' risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`, then
#' it further includes `n_factors`-dimensional vector of factor risk
#' premia standard errors in `"standard_errors"`.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute KRS factor risk premia and their standard errors
#' frp = FRP(returns, factors, include_standard_errors = TRUE)
#'
#' @export
FRP = function(
  returns,
  factors,
  misspecification_robust = TRUE,
  include_standard_errors = FALSE,
  check_arguments = TRUE
) {

  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`misspecification_robust` must be boolean" = is.logical(misspecification_robust))
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))

  }

  return(.Call(`_intrinsicFRP_FRPCpp`,
    returns,
    factors,
    misspecification_robust,
    include_standard_errors
  ))

}
