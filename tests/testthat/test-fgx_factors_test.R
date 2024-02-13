test_that("Test FGXFactorsTest", {

  factors = factors
  returns = returns[,-1]

  # Function returns a list
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7])
  expect_type(output, "list")

  # Function returns a list with three elements
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7])
  expect_length(output, 3)

  # SDF coefficients are numeric
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7])
  expect_type(output$sdf_coefficients, "double")

  # Standard errors are numeric
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7])
  expect_type(output$standard_errors, "double")

  # Number of SDF coefficients equals number of new factors
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7])
  expect_length(output$sdf_coefficients, ncol(factors[, 5:7]))

  # Number of standard errors equals number of new factors
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7])
  expect_length(output$standard_errors, ncol(factors[, 5:7]))

  # Error is thrown for missing values in new_factors
  expect_error(FGXFactorsTest(returns, factors[, 2:4], cbind(factors[, 5:7], NA)))

  # Error is thrown for non-numeric new_factors
  expect_error(FGXFactorsTest(returns, factors[, 2:4], as.character(factors[, 5:7])))

  # Error is thrown for new_factors with more variables than observations
  expect_error(FGXFactorsTest(returns, factors[, 2:4], factors[1:10, 5:7]))

  # Error is thrown for non-numeric n_folds
  expect_error(FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7], nfolds = "a"))

  # Function runs without error when check_arguments = FALSE
  output <- FGXFactorsTest(returns, factors[, 2:4], factors[, 5:7], check_arguments = FALSE)
  expect_type(output, "list")

  # Error is thrown when check_arguments = FALSE and new_factors contain missing values
  expect_error(FGXFactorsTest(returns, factors[, 2:4], cbind(factors[, 5:7], NA), check_arguments = FALSE))

})

