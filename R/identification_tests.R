# Author: Alberto Quaini

#####################################
######  Identification tests ########
#####################################

#' @title Asset Pricing Model Identification via Chen-Fang (2019) Beta Rank Test
#'
#' @name ChenFang2019BetaRankTest
#' @description Tests the null hypothesis of reduced rank in the matrix of regression
#' loadings for test asset excess returns on risk factors using the Chen-Fang (2019)
#' <doi:10.3982/QE1139>
#' beta rank test. The test applies the Kleibergen-Paap (2006) <doi:10.1016/j.jeconom.2005.02.011>
#'  iterative rank test
#' for initial rank estimation when `target_level_kp2006_rank_test > 0`, with an
#' adjustment to `level = target_level_kp2006_rank_test / n_factors`. When
#' `target_level_kp2006_rank_test <= 0`, the number of singular values above
#' `n_observations^(-1/4)` is used instead. It presumes that the number of factors
#' is less than the number of returns (`n_factors < n_returns`).
#' All the details can be found in Chen-Fang (2019)
#' <doi:10.3982/QE1139>.
#'
#' @param returns Matrix of test asset excess returns with dimensions `n_observations x n_returns`.
#' @param factors Matrix of risk factors with dimensions `n_observations x n_factors`.
#' @param n_bootstrap The number of bootstrap samples to use in the Chen-Fang (2019) test.
#' Defaults to 500 if not specified.
#' @param target_level_kp2006_rank_test The significance level for the Kleibergen-Paap (2006)
#' rank test used for initial rank estimation. If set above 0, it indicates the level for this
#' estimation within the Chen-Fang (2019) rank test. If set at 0 or negative, the initial rank
#' estimator defaults to the count of singular values exceeding `n_observations^(-1/4)`.
#' The default value is `0.05` to account for multiple testing.
#' @param check_arguments Logical flag to determine if input arguments should be checked for validity.
#' Default is `TRUE`.
#'
#' @return A list containing the Chen-Fang (2019) rank statistic and the associated p-value.
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
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

  # check function arguments
  if (check_arguments) {
    CheckData(returns, factors)
    stopifnot("`n_bootstrap` must be numeric" = is.numeric(n_bootstrap))
    stopifnot("`n_bootstrap` must be greater than 0" = n_bootstrap > 0)
    stopifnot("`target_level_kp2006_rank_test` must be numeric" = is.numeric(target_level_kp2006_rank_test))
    stopifnot("`target_level_kp2006_rank_test` must be between 0 and 1" = (target_level_kp2006_rank_test >= 0.) & (target_level_kp2006_rank_test <= 1.))
  }

  # compute the Chen-Fang (2019) test statistic and p-value
  return(.Call(`_intrinsicFRP_ChenFang2019BetaRankTestCpp`,
    returns,
    factors,
    n_bootstrap,
    target_level_kp2006_rank_test
  ))

}

#' @title Asset Pricing Model Identification via Iterative Kleibergen-Paap 2006 Beta Rank Test
#'
#' @name IterativeKleibergenPaap2006BetaRankTest
#' @description Evaluates the rank of regression loadings in an asset pricing model using the
#' iterative Kleibergen-Paap (2006) <doi:10.1016/j.jeconom.2005.02.011> beta rank test.
#' It systematically tests the null hypothesis
#' for each potential rank `q = 0, ..., n_factors - 1` and estimates the rank as the smallest `q`
#' that has a p-value below the significance level, adjusted for the number of factors.
#' The function presupposes more returns than factors (`n_factors < n_returns`).
#' All the details can be found in Kleibergen-Paap (2006) <doi:10.1016/j.jeconom.2005.02.011>.
#'
#' @param returns A matrix of test asset excess returns with dimensions `n_observations x n_returns`.
#' @param factors A matrix of risk factors with dimensions `n_observations x n_factors`.
#' @param target_level A numeric value specifying the significance level for the test. For each
#' hypothesis test `H: rank(beta) = q`, the significance level is adjusted to
#' `target_level / n_factors`. The default is `0.05`.
#' @param check_arguments Logical flag indicating whether to perform internal checks of the
#' function's arguments. Defaults to `TRUE`.
#'
#' @return A list containing estimates of the regression loading rank and the associated
#' iterative Kleibergen-Paap 2006 beta rank statistics and p-values for each `q`.
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
IterativeKleibergenPaap2006BetaRankTest = function(
  returns,
  factors,
  target_level = 0.05,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {
    CheckData(returns, factors)
    stopifnot("`target_level` must be numeric" = is.numeric(target_level))
    stopifnot("`target_level` must be between 0 and 1" = (target_level >= 0.) & (target_level <= 1.))
  }

  # compute the Kleibergen-Paap (2006) test
  return(.Call(`_intrinsicFRP_IterativeKleibergenPaap2006BetaRankTestCpp`,
    returns,
    factors,
    target_level
  ))

}
