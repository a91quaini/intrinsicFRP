test_that("Test GKRFactorScreening", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  target_level <- 0.05
  hac_prewhite <- FALSE

  # Test basic functionality without errors
  expect_no_error(GKRFactorScreening(returns, factors, target_level, hac_prewhite))
  expect_no_error(GKRFactorScreening(
    returns,
    matrix(stats::rnorm(nrow(factors) * n_factors), nrow(factors), n_factors),
    target_level,
    hac_prewhite
  ))

  # Test with prewhitening
  expect_no_error(GKRFactorScreening(returns, factors, target_level, TRUE))

  # Test error handling for incorrect dimensions
  expect_error(GKRFactorScreening(t(returns), factors, target_level, hac_prewhite))
  expect_error(GKRFactorScreening(returns, t(factors), target_level, hac_prewhite))

  # Test errors for wrong input types
  expect_error(GKRFactorScreening(c(), factors, target_level, hac_prewhite))
  expect_error(GKRFactorScreening(returns, c(), target_level, hac_prewhite))

  # Test for invalid target levels
  expect_error(GKRFactorScreening(returns, factors, -0.1, hac_prewhite))
  expect_error(GKRFactorScreening(returns, factors, 1.1, hac_prewhite))

  # More comprehensive tests can include checking the structure of the output
  screening_result <- GKRFactorScreening(returns, factors, target_level, hac_prewhite)
  expect_type(screening_result, "list")
  expect_true(all(c("sdf_coefficients", "standard_errors", "t_statistics", "selected_factor_indices") %in% names(screening_result)))

  # Validate the length of the output lists
  expect_no_error(length(screening_result$sdf_coefficient) <= n_factors)
  expect_no_error(length(screening_result$standard_errors) <= n_factors)
  expect_no_error(length(screening_result$squared_t_stat) <= n_factors)
  expect_no_error(length(screening_result$selected_factor_indices) <= n_factors)

})
