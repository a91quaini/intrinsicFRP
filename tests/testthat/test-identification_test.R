test_that("Test ChenFang2019BetaRankTest", {

  factors = factors[,-1]
  returns = returns[,-1]


  useless_factors = matrix(rnorm(nrow(factors) * 10), nrow(factors), 10)

  expect_no_error(ChenFang2019BetaRankTest(returns, factors))
  expect_no_error(ChenFang2019BetaRankTest(
    returns,
    factors,
    n_bootstrap = 400,
    target_level_kp2006_rank_test = 0.05
  ))
  expect_no_error(ChenFang2019BetaRankTest(
    returns,
    cbind(factors, useless_factors),
    n_bootstrap = 400,
    target_level_kp2006_rank_test = 0.99
  ))
  expect_no_error(ChenFang2019BetaRankTest(
    returns,
    cbind(factors, useless_factors),
    n_bootstrap = 400,
    target_level_kp2006_rank_test = 1.e-4
  ))
  expect_error(
    ChenFang2019BetaRankTest(returns, factors, n_bootstrap = "r")
  )
  expect_error(
    ChenFang2019BetaRankTest(returns, factors, target_level_kp2006_rank_test = "r")
  )

})
