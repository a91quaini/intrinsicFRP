# Author: Alberto Quaini

#############################################
######  HJ misspecification distance ########
#############################################


#' @title Compute the HJ asset pricing model misspecification distance.
#'
#' @name HJMisspecificationDistance
#' @description Computes the Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>
#' squared model misspecification distance:
#' `square_distance = min_{d} (E[R] - Cov[R,F] * d)' * V[R]^{-1} * (E[R] - Cov[R,F] * d)`,
#' where `R` denotes test asset excess returns and `F` risk factors,
#' and computes the associated confidence interval.
#' This model misspecification distance is a modification of the prominent
#' Hansen-Jagannathan (1997) <doi:10.1111/j.1540-6261.1997.tb04813.x>
#' distance, adapted to the use of excess returns for the test asset, and a
#' SDF that is a linear function of demeaned factors.
#' Clearly, computation of the confidence interval is obtained by means of an
#' asymptotic analysis under potentially misspecified models, i.e.,
#' without assuming correct model specification.
#' Details can be found in Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>.
#'
#' @param returns A `n_observations x n_returns` matrix of test asset excess returns.
#' @param factors A `n_observations x n_factors` matrix of risk factors.
#' @param ci_coverage A number indicating the confidence interval coverage
#' probability. Default is `0.95`.
#' @param hac_prewhite A boolean indicating if the series needs pre-whitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param check_arguments A boolean: `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return @return A list containing the squared misspecification-robust HJ
#' distance in `squared_distance`, and the lower and upper confidence bounds
#' in `lower_bound` and `upper_bound`, respectively.
#'
#' @examples
#' # Import package data on 6 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # Compute the HJ model misspecification distance
#' hj_test = HJMisspecificationDistance(returns, factors)
#'
#' @export
HJMisspecificationDistance = function(
  returns,
  factors,
  ci_coverage = 0.95,
  hac_prewhite = FALSE,
  check_arguments = TRUE
) {

  # Check the function arguments if check_arguments is TRUE
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`ci_coverage` must be numeric" = is.numeric(ci_coverage))
    stopifnot("`ci_coverage` must be between  0 and 1" = (ci_coverage >= 0.) & (ci_coverage <= 1.))
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))

  }

  # Call the C++ implementation of the HJ Misspecification distance
  return(.Call(`_intrinsicFRP_HJMisspecificationDistanceCpp`,
    returns,
    factors,
    ci_coverage,
    hac_prewhite
  ))

}
