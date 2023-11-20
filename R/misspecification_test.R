# Author: Alberto Quaini

#########################################
######  HJ misspecification test ########
#########################################


#' @title Asset Prciing HJ model misspecification test
#'
#' @description Computes the Hansen-Jagannathan (1997) <doi:10.1111/j.1540-6261.1997.tb04813.x>
#' model misspecification statistic and
#' p-value of an asset pricing model from test asset excess returns `R` and
#' risk factors `F`. The statistic is defined as:
#' `HJDISTANCE = min_{d} (E[R] - beta * d)' * Var[R]^{-1} * (E[R] - beta * d)`
#' where `beta = Cov[R, F] * Var[F]^{-1}` are the regression coefficients of
#' test asset excess returns `R` on risk factors `F`.
#' Detailed computations and p-value calculations can be found in
#' Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>.
#'
#' @param returns An `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors An `n_observations x n_factors`-dimensional matrix of risk
#' factors.
#' @param hac_prewhite A boolean indicating if the series needs prewhitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param check_arguments A boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A list containing the squared standardized HJ test statistic and
#' the corresponding p-value.
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
  hac_prewhite = FALSE,
  check_arguments = TRUE
) {

  # Check the function arguments if check_arguments is TRUE
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))

  }

  # Call the C++ implementation of the HJ Misspecification Test
  return(.Call(`_intrinsicFRP_HJMisspecificationTestCpp`,
    returns,
    factors,
    hac_prewhite
  ))

}
