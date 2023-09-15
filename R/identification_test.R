# Author: Alberto Quaini

####################################
######  Identification test ########
####################################

#' Compute the Chen Fang 2019 beta rank test
#'
#' @name ChenFang2019BetaRankTest
#' @description Computes the Chen fang 2019 rank statistic and p-value of the
#' null that the matrix of regression loadings of test asset excess returns on
#' risk factors has reduced rank. If `target_level_kp2006_rank_test > 0`,
#' it uses the iterative Kleibergen Paap 2006 rank test to estimate the initial
#' rank, with `level = target_level_kp2006_rank_test / n_factors`. If
#' `target_level_kp2006_rank_test <= 0`, the initial rank estimator is taken to
#' be the number of singular values above `n_observations^(-1/4)`. It assumes
#' n_factors < n_returns.
#'
#' @param returns `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors `n_observations x n_factors`-dimensional matrix of risk
#' factors.
#' @param n_bootstrap numeric integer indicating the number of bootstrap
#' samples used to compute the Chen fang 2019 test. Default is `500`.
#' @param target_level_kp2006_rank_test numeric level of the Kleibergen Paap
#' 2006 rank test. If it is strictly grater than zero, then the iterative
#' Kleibergen Paap 2006 rank test at `level = target_level_kp2006_rank_test /
#' n_factors` is used to compute an initial estimator of the rank of the factor
#' loadings in the Chen Fang 2019 rank test. Otherwise, the initial rank
#' estimator is taken to be the number of singular values above
#' `n_observations^(-1/4)`. Default is `0.05` (as correction for multiple
#' testing).
#' @param check_arguments boolean `TRUE` if you want to check function arguments;
#' `FALSE` otherwise. Default is `TRUE`.
#'
#' @return a list containing the Chen fang 2019 rank statistic and the
#' corresponding p-value.
#'
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute the model identification test
#' hj_test = ChenFang2019BetaRankTest(returns, factors)
#'
#' @export
ChenFang2019BetaRankTest = function(
  returns,
  factors,
  n_bootstrap = 500,
  target_level_kp2006_rank_test = 0.05,
  check_arguments = TRUE
) {

  if (check_arguments) {
    CheckData(returns, factors)
    stopifnot("`n_bootstrap` must be numeric" = is.numeric(n_bootstrap))
    stopifnot("`n_bootstrap` must be greater than 0" = n_bootstrap > 0)
    stopifnot("`target_level_kp2006_rank_test` must be numeric" = is.numeric(target_level_kp2006_rank_test))
    stopifnot("`target_level_kp2006_rank_test` must be between 0 and 1" = (target_level_kp2006_rank_test >= 0.) & (target_level_kp2006_rank_test <= 1.))
  }

  return(.Call(`_intrinsicFRP_ChenFang2019BetaRankTestCpp`,
    returns,
    factors,
    n_bootstrap,
    target_level_kp2006_rank_test
  ))

}
