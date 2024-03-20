# Author: Alberto Quaini

###################################
######  SDF coefficients ##########
###################################

#' @title SDF Coefficients
#'
#' @name SDFCoefficients
#' @description Computes the misspecification-robust SDF coefficients of
#' Gospodinov-Kan-Robotti (2014) <https:#'doi.org/10.1093/rfs/hht135>:
#' `GKRSDFcoefficients = (C' * V[R]^{-1} * C)^{-1} * C' * V[R]^{-1} * E[R]`
#' from data on factors `F` and test
#' asset excess returns `R`.
#' These notions of SDF coefficients minimize pricing errors:
#' `argmin_{d} (E[R] - Cov[R,F] * d)' * V[R]^{-1} * (E[R] - Cov[R,F] * d)`.
#' Optionally computes the corresponding
#' heteroskedasticity and autocorrelation robust standard errors using the
#' Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the
#' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#'
#' @param returns A `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors A `n_observations x n_factors`-dimensional matrix of factors.
#' @param include_standard_errors A boolean: `TRUE` if you want to compute the
#' SDF coefficients' HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param hac_prewhite A boolean indicating if the series needs pre-whitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param check_arguments A boolean: `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A list containing `n_factors`-dimensional vector of SDF coefficients
#' in `"sdf_coefficients"`; if `include_standard_errors = TRUE`, then
#' it further includes `n_factors`-dimensional vector of SDF coefficients'
#' standard errors in `"standard_errors"`;
#'
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute GKR SDF coefficients and their standard errors
#' frp = SDFCoefficients(returns, factors, include_standard_errors = TRUE)
#'
#' @export
SDFCoefficients = function(
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

  # compute the FRP estimate and, eventually, their standard errors
  return(.Call(`_intrinsicFRP_SDFCoefficientsCpp`,
    returns,
    factors,
    include_standard_errors,
    hac_prewhite
  ))

}
