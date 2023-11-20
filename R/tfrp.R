# Author: Alberto Quaini

#######################################
######  Tradable Factor Risk Premia ###
#######################################

#' Compute tradable factor risk premia from data
#'
#' @name TFRP
#' @description Computes tradable factor risk premia from data on factors `F` and
#' test asset excess returns `R`:
#' `TFRP = Cov[F, R] * Var[R]^{-1} * E[R]`.
#' Optionally computes the corresponding heteroskedasticity and autocorrelation
#' robust standard errors using the Newey-West (1994) <doi:10.2307/2297912>
#' plug-in procedure to select the number of relevant lags, i.e.,
#' `n_lags = 4 * (n_observations/100)^(2/9)`.
#' All details are found in Quaini-Trojani-Yuan (2023) <doi:10.2139/ssrn.4574683>.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of factors.
#' @param include_standard_errors boolean `TRUE` if you want to compute the
#' tradable factor risk premia HAC standard errors; `FALSE` otherwise. Default
#' is `FALSE`.
#' @param hac_prewhite A boolean indicating if the series needs prewhitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param check_arguments boolean `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing `n_factors`-dimensional vector of tradable factor
#' risk premia in `"risk_premia"`; if `include_standard_errors=TRUE`, then
#' it further includes `n_factors`-dimensional vector of tradable factor risk
#' premia standard errors in `"standard_errors"`.
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute tradable factor risk premia and their standard errors
#' tfrp = TFRP(returns, factors, include_standard_errors = TRUE)
#'
#' @export
TFRP = function(
  returns,
  factors,
  include_standard_errors = FALSE,
  hac_prewhite = FALSE,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))

  }

  # compute the TFRP estimate and, eventually, their standard errors
  return(.Call(`_intrinsicFRP_TFRPCpp`,
    returns,
    factors,
    include_standard_errors,
    hac_prewhite
  ))

}

