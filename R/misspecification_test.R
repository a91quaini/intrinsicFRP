# Author: Alberto Quaini

#########################################
######  HJ misspecification test ########
#########################################


#' Compute factor risk premia
#'
#' @name HJMisspecificationTest
#' @description Computes the Hansen-Jagannatan misspecification statistic and
#' p-value of an asset pricing model from test asset excess returns and
#' risk factors.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of risk
#' factors.
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing the HJ test statistic and the corresponding
#' p-value.
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute the HJ model misspecification test
#' hj_test = HJMisspecificationTest(returns, factors)
#'
#' @export
HJMisspecificationTest = function(returns, factors, check_arguments = TRUE) {

  if (check_arguments) {CheckData(returns, factors)}

  return(.Call(`_intrinsicFRP_HJMisspecificationTestCpp`,
    returns,
    factors
  ))

}
