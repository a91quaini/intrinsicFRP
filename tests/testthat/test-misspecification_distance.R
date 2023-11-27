test_that("Test HJMisspecificationDistance", {

  factors = factors[,-1]
  returns = returns[,-1]

  # Test if the function works without errors for valid inputs.
  expect_no_error(HJMisspecificationDistance(returns, factors))

  # Test if prewhite works
  expect_no_error(HJMisspecificationDistance(returns, factors, hac_prewhite = TRUE))

  # Test if the function correctly throws an error when an empty matrix is passed as 'returns'.
  expect_error(HJMisspecificationDistance(
    matrix(numeric(0), nrow = nrow(returns), ncol = ncol(returns)),
    factors
  ))

  # Test if the function correctly throws an error when an empty matrix is passed as 'factors'.
  expect_error(HJMisspecificationDistance(
    returns,
    matrix(numeric(0), nrow = nrow(factors), ncol = ncol(factors))
  ))

  # Test if the function correctly throws an error when 'returns' has fewer rows than 'factors'.
  expect_error(HJMisspecificationDistance(
    returns[1:(nrow(returns)-5),],
    factors
  ))

  # Test if the function correctly throws an error when 'factors' has fewer rows than 'returns'.
  expect_error(HJMisspecificationDistance(
    returns,
    factors[1:(nrow(factors)-5),]
  ))

  # Test if the function is consistent in its output given the same inputs.
  expect_equal(
    HJMisspecificationDistance(returns, factors),
    HJMisspecificationDistance(returns, factors)
  )

  # Test if function correctly throws an error when 'hac_prewhite' is not boolean.
  expect_error(HJMisspecificationDistance(returns, factors, hac_prewhite = "c"))

  # Test different values of ci_coverage within the valid range
  expect_no_error(HJMisspecificationDistance(returns, factors, ci_coverage = 0.90))
  expect_no_error(HJMisspecificationDistance(returns, factors, ci_coverage = 0.99))

  # Test ci_coverage outside the valid range
  expect_error(HJMisspecificationDistance(returns, factors, ci_coverage = -0.1))
  expect_error(HJMisspecificationDistance(returns, factors, ci_coverage = 1.1))

  # Test for edge cases with minimal dimensions
  expect_no_error(HJMisspecificationDistance(matrix(rnorm(10), ncol = 1), matrix(rnorm(10), ncol = 1)))

  # Test for non-numeric matrices
  expect_error(HJMisspecificationDistance(matrix("a", nrow = 10, ncol = 2), factors))
  expect_error(HJMisspecificationDistance(returns, matrix("b", nrow = 10, ncol = 2)))

  # Test the structure of the output
  result <- HJMisspecificationDistance(returns, factors)
  expect_true(is.list(result))
  expect_true(all(c("squared_distance", "lower_bound", "upper_bound") %in% names(result)))

  # Test for handling of NA, NaN, Inf values
  returns_with_na <- returns
  returns_with_na[1,1] <- NA
  expect_error(HJMisspecificationDistance(returns_with_na, factors))

  # Stress test with large dataset
  large_returns <- matrix(rnorm(10000), ncol = 10)
  large_factors <- matrix(rnorm(10000), ncol = 10)
  expect_no_error(HJMisspecificationDistance(large_returns, large_factors))

  # Check for side effects on inputs
  original_returns <- returns
  HJMisspecificationDistance(returns, factors)
  expect_identical(returns, original_returns)


})

