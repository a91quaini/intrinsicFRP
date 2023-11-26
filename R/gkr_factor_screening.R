# Author: Alberto Quaini

###################################
######  GKRFactorScreeening #######
###################################

#' @title Factor screening procedure of Gospodinov-Kan-Robotti (2014)
#'
#' @name GKRFactorScreening
#' @description Performs the factor screening procedure of
#' Gospodinov-Kan-Robotti (2014) <doi:10.2139/ssrn.2579821>, which is
#' an iterative model screening procedure
#' based on the sequential removal of factors associated with the smallest insignificant
#' t-test of a nonzero misspecification-robust SDF coefficient. The significance threshold for the
#' absolute t-test is set to `target_level / n_factors`,
#' where n_factors indicates the number of factors in the model at the current iteration;
#' that is, it takes care of the multiple testing problem via a conservative
#' Bonferroni correction. Standard errors are computed with the
#' heteroskedasticity and autocorrelation using the Newey-West (1994)
#' <doi:10.2307/2297912> estimator, where the number of lags
#' is selected using the Newey-West plug-in procedure:
#' `n_lags = 4 * (n_observations/100)^(2/9)`.
#' For the standard error computations, the function allows to internally
#' pre-whiten the series by fitting a VAR(1),
#' i.e., a vector autoregressive model of order 1.
#' All the details can be found in Gospodinov-Kan-Robotti (2014) <doi:10.2139/ssrn.2579821>.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of risk
#' factors.
#' @param target_level Number specifying the target significance threshold for the
#' tests underlying the GKR factor screening procedure.
#' To account for the multiple testing problem, the significance threshold for the
#' absolute t-test is given by `target_level_gkr2014_screening / n_factors`,
#' where n_factors indicate the number of factors in the model at the current iteration.
#' Default is `0.05`.
#' @param hac_prewhite A boolean indicating if the series needs prewhitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param check_arguments boolean `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A list contaning the selected GKR SDF coefficients in `SDF_coefficients`,
#' their standard errors in `standard_errors`,
#' t-statistics in `t_statistics` and indices in the columns of the factor matrix `factors`
#' supplied by the user in `selected_factor_indices`.
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # Perform the GKR factor screening procedure
#' screen = GKRFactorScreening(returns, factors)
#'
#' @export
GKRFactorScreening = function(
  returns,
  factors,
  target_level = 0.05,
  hac_prewhite = FALSE,
  check_arguments = TRUE
) {

  # Check function arguments.
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`target_level` must be numeric" = is.numeric(target_level))
    stopifnot("`target_level` must be between 0 and 1" = (target_level >= 0.) & (target_level <= 1.))
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))

  }

  # Perform the GKR factor screening procedure.
  results = .Call(`_intrinsicFRP_GKRFactorScreeningCpp`,
    returns,
    factors,
    target_level,
    hac_prewhite
  )

  # Transform c++ indices (starting at 0) into R indices (starting at 1).
  results$selected_factor_indices = results$selected_factor_indices + 1

  return(results)

}
