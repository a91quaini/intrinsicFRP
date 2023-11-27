test_that("Test ChenFang2019BetaRankTest and
  IterativeKleibergenPaap2006BetaRankTest", {

  factors = factors[,-1]
  returns = returns[,-1]

  # Creating a matrix of random noise to use as useless factors.
  useless_factors = matrix(rnorm(nrow(factors) * 10), nrow(factors), 10)

  ## ChenFang2019BetaRankTest

  # Testing the basic functionality of ChenFang2019BetaRankTest without errors.
  expect_no_error(ChenFang2019BetaRankTest(returns, factors))

  # Testing ChenFang2019BetaRankTest with additional bootstrap and target level arguments.
  expect_no_error(ChenFang2019BetaRankTest(
    returns,
    factors,
    n_bootstrap = 400,
    target_level_kp2006_rank_test = 0.05
  ))

  # Testing the robustness of ChenFang2019BetaRankTest with added noise factors.
  expect_no_error(ChenFang2019BetaRankTest(
    returns,
    cbind(factors, useless_factors),
    n_bootstrap = 400,
    target_level_kp2006_rank_test = 0.99
  ))

  # Testing ChenFang2019BetaRankTest with a small target level.
  expect_no_error(ChenFang2019BetaRankTest(
    returns,
    cbind(factors, useless_factors),
    n_bootstrap = 400,
    target_level_kp2006_rank_test = 1.e-4
  ))

  # Testing error handling in ChenFang2019BetaRankTest for invalid n_bootstrap and target_level arguments.
  expect_error(
    ChenFang2019BetaRankTest(returns, factors, n_bootstrap = "r")
  )
  expect_error(
    ChenFang2019BetaRankTest(returns, factors, target_level_kp2006_rank_test = "r")
  )

  ## IterativeKleibergenPaap2006BetaRankTest

  # Testing the basic functionality of IterativeKleibergenPaap2006BetaRankTest without errors.
  expect_no_error(IterativeKleibergenPaap2006BetaRankTest(returns, factors))

  # Testing IterativeKleibergenPaap2006BetaRankTest with a specified target level.
  expect_no_error(IterativeKleibergenPaap2006BetaRankTest(
    returns,
    factors,
    target_level = 0.05
  ))

  # Testing the robustness of IterativeKleibergenPaap2006BetaRankTest with added noise factors.
  expect_no_error(IterativeKleibergenPaap2006BetaRankTest(
    returns,
    cbind(factors, useless_factors),
    target_level = 0.99
  ))

  # Testing IterativeKleibergenPaap2006BetaRankTest with a small target level.
  expect_no_error(IterativeKleibergenPaap2006BetaRankTest(
    returns,
    cbind(factors, useless_factors),
    target_level = 1.e-4
  ))

  # Testing error handling in IterativeKleibergenPaap2006BetaRankTest for invalid target_level argument.
  expect_error(
    IterativeKleibergenPaap2006BetaRankTest(returns, factors, target_level = "r")
  )

  # Test with random factors
  random_factors = matrix(rnorm(nrow(returns) * (ncol(returns) - 1)), ncol = ncol(returns) - 1)
  expect_no_error(ChenFang2019BetaRankTest(returns, random_factors))

  # Test output structure and type
  cf_result <- ChenFang2019BetaRankTest(returns, factors)
  expect_true(is.list(cf_result))
  expect_true(all(c("statistic", "p-value") %in% names(cf_result)))

  # Test edge cases for bootstrap and target level
  expect_no_error(ChenFang2019BetaRankTest(returns, factors, n_bootstrap = 1))
  expect_no_error(ChenFang2019BetaRankTest(returns, factors, target_level_kp2006_rank_test = 0))
  expect_no_error(ChenFang2019BetaRankTest(returns, factors, target_level_kp2006_rank_test = 1))

  # Test for consistency of output
  expect_equal(
    ChenFang2019BetaRankTest(returns, factors),
    ChenFang2019BetaRankTest(returns, factors)
  )

  # Test for handling of NA or Inf values
  returns_with_na <- returns
  returns_with_na[1,1] <- NA
  expect_error(ChenFang2019BetaRankTest(returns_with_na, factors))

  # Test for input matrices with one column
  expect_no_error(ChenFang2019BetaRankTest(returns[,1:2,drop=FALSE], factors[,1,drop=FALSE]))

  # Test for incorrect dimensions
  expect_error(ChenFang2019BetaRankTest(returns, factors[-(1:5),]))

  # Test for invalid data types
  expect_error(ChenFang2019BetaRankTest(returns, as.character(factors)))
  })
