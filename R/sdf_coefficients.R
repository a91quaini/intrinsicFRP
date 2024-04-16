# Author: Alberto Quaini

###################################
######  SDF coefficients ##########
###################################

#' @title SDF Coefficients
#'
#' @name SDFCoefficients
#' @description Computes the SDF coefficients of Fama-MachBeth (1973) <doi:10.1086/260061>
#' `FMSDFcoefficients = (C' * C)^{-1} * C' * E[R]`
#' or the misspecification-robust SDF coefficients of
#' Gospodinov-Kan-Robotti (2014) <doi:10.1093/rfs/hht135>:
#' `GKRSDFcoefficients = (C' * V[R]^{-1} * C)^{-1} * C' * V[R]^{-1} * E[R]`
#' from data on factors `F` and test asset excess returns `R`.
#' These notions of SDF coefficients minimize pricing errors:
#' `argmin_{d} (E[R] - Cov[R,F] * d)' * W * (E[R] - Cov[R,F] * d)`,
#' with `W=I`, i.e., the identity, and `W=V[R]^{-1}`, respectively.
#' Optionally computes the corresponding
#' heteroskedasticity and autocorrelation robust standard errors (accounting
#' for a potential model misspecification) using the
#' Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the
#' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#'
#' @param returns A `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors A `n_observations x n_factors`-dimensional matrix of factors.
#' @param misspecification_robust A boolean: `TRUE` for the
#' "misspecification-robust" Kan-Robotti-Shanken (2013) GLS approach using the
#' inverse covariance matrix of returns; `FALSE` for standard Fama-MacBeth
#' risk premia. Default is `TRUE`.
#' @param include_standard_errors A boolean: `TRUE` if you want to compute the
#' SDF coefficients' HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param hac_prewhite A boolean indicating if the series needs pre-whitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param target_level_gkr2014_screening A number indicating the target level of
#' the tests underlying the factor screening procedure in Gospodinov-Kan-Robotti
#' (2014). If it is zero, then no factor screening procedure is
#' implemented. Otherwise, it implements an iterative screening procedure
#' based on the sequential removal of factors associated with the smallest insignificant
#' t-test of a nonzero SDF coefficient. The threshold for the absolute t-test is
#' `target_level_gkr2014_screening / n_factors`, where n_factors indicate the
#' number of factors in the model at the current iteration. Default is `0.`, i.e.,
#' no factor screening.
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
  misspecification_robust = TRUE,
  include_standard_errors = FALSE,
  hac_prewhite = FALSE,
  target_level_gkr2014_screening = 0.,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`misspecification_robust` must be boolean" = is.logical(misspecification_robust))
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))
    stopifnot("`target_level_gkr2014_screening` must be numeric" = is.numeric(target_level_gkr2014_screening))
    stopifnot("`target_level_gkr2014_screening` must be between 0 and 1" = (target_level_gkr2014_screening >= 0.) & (target_level_gkr2014_screening <= 1.))

  }

  # compute the FRP estimate and, eventually, their standard errors
  output = .Call(`_intrinsicFRP_SDFCoefficientsCpp`,
    returns,
    factors,
    misspecification_robust,
    include_standard_errors,
    hac_prewhite,
    target_level_gkr2014_screening
  )

  if (target_level_gkr2014_screening > 0) {
    # Transform c++ indices (starting at 0) into R indices (starting at 1).
    output$selected_factor_indices = output$selected_factor_indices + 1
  }

  return(output)

}
