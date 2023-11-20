test_that("Test HJMisspecificationTest", {

  factors = factors[,-1]
  returns = returns[,-1]

  # Test if the function works without errors for valid inputs.
  expect_no_error(HJMisspecificationTest(returns, factors))

  # Test if prewhite works
  expect_no_error(HJMisspecificationTest(returns, factors, hac_prewhite = TRUE))

  # Test if the function correctly throws an error when an empty matrix is passed as 'returns'.
  expect_error(HJMisspecificationTest(
    matrix(numeric(0), nrow = nrow(returns), ncol = ncol(returns)),
    factors
  ))

  # Test if the function correctly throws an error when an empty matrix is passed as 'factors'.
  expect_error(HJMisspecificationTest(
    returns,
    matrix(numeric(0), nrow = nrow(factors), ncol = ncol(factors))
  ))

  # Test if the function correctly throws an error when 'returns' has fewer rows than 'factors'.
  expect_error(HJMisspecificationTest(
    returns[1:(nrow(returns)-5),],
    factors
  ))

  # Test if the function correctly throws an error when 'factors' has fewer rows than 'returns'.
  expect_error(HJMisspecificationTest(
    returns,
    factors[1:(nrow(factors)-5),]
  ))

  # Test if the function is consistent in its output given the same inputs.
  expect_equal(
    HJMisspecificationTest(returns, factors),
    HJMisspecificationTest(returns, factors)
  )

  # Test if function correctly throws an error when 'hac_prewhite' is not boolean.
  expect_error(HJMisspecificationTest(returns, factors, hac_prewhite = "c"))


})
