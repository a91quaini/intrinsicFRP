# Author: Alberto Quaini

#########################################
######  HJ misspecification test ########
#########################################


#' @title Compute the HJ asset pricing model misspecification test.
#'
#' @name HJMisspecificationTest
#' @description Computes and tests the Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>
#' squared model misspecification distance:
#' `sqhj_distance = min_{d} (E[R] - Cov[R,F] * d)' * V[R]^{-1} * (E[R] - Cov[R,F] * d)`,
#' where `R` denotes test asset excess returns and `F` risk factors.
#' This model misspecification distance is a modification of the prominent
#' Hansen-Jagannathan (1997) <doi:10.1111/j.1540-6261.1997.tb04813.x>
#' distance, adapted to the use of excess returns for the test asset, and a
#' SDF that is a linear function of demeaned factors.
#' The Null Hypothesis of the test is:
#' `H0: sqhj_distance = sqhj_distance_null_value`,
#' where `sqhj_distance_null_value` is the user-supplied hypothesized value of
#' the squared model misspecification distance.
#' Computation of the p-values is obtained with asymptotic analysis under
#' misspecified models; that is, without assuming correct specification.
#' Details can be found in Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>.
#'
#' @param returns A `n_observations x n_returns` matrix of test asset excess returns.
#' @param factors A `n_observations x n_factors` matrix of risk factors.
#' @param sqhj_distance_null_value A number indicating the testing value under
#' the null hypothesisfor the squared HJ model-misspecification distance.
#' Default is zero.
#' @param hac_prewhite A boolean indicating if the series needs pre-whitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param check_arguments A boolean: `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return @return A list containing the squared misspecification-robust HJ
#' distance, the associated standardized test statistic
#' `T * (sqhj_distance - sqhj_distance_null_value)^2 / AsySdErr(sqhj_distance)`
#' and the corresponding p-value.
#'
#' @examples
#' # Import package data on 6 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # Compute the HJ model misspecification test
#' hj_test = HJMisspecificationTest(returns, factors)
#'
#' @export
HJMisspecificationTest = function(
  returns,
  factors,
  sqhj_distance_null_value = 0.,
  hac_prewhite = FALSE,
  check_arguments = TRUE
) {

  # Check the function arguments if check_arguments is TRUE
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`sqhj_distance_null_value` must be numeric" = is.numeric(sqhj_distance_null_value))
    stopifnot("`sqhj_distance_null_value` must be greater or equal to zero" = sqhj_distance_null_value >= 0.)
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))

  }

  # Call the C++ implementation of the HJ Misspecification Test
  return(.Call(`_intrinsicFRP_HJMisspecificationTestCpp`,
    returns,
    factors,
    sqhj_distance_null_value,
    hac_prewhite
  ))

}
