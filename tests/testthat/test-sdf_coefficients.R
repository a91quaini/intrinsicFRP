test_that("Test SDFCoefficients", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  # Calculating factor loadings using Fama-MacBeth procedure.
  covariance_returns_factors = stats::cov(returns, factors)
  variance_returns = stats::cov(returns)
  mean_returns = matrix(colMeans(returns), n_returns, 1)

  # Calculating SDF coefficients using the Gospodinov-Kan-Robotti approach.
  cov_fac_ret_var_ret_inv = t(solve(
    variance_returns,
    covariance_returns_factors
  ))
  sdf_coeffs = solve(
    cov_fac_ret_var_ret_inv %*% covariance_returns_factors,
    cov_fac_ret_var_ret_inv
  ) %*% mean_returns

  # Testing basic functionality of SDFCoefficients without errors.
  expect_no_error(SDFCoefficients(returns, factors))

  # Testing SDFCoefficients with standard errors included.
  expect_no_error(SDFCoefficients(returns, factors, include_standard_errors = TRUE))

  # Test if prewhite works
  expect_no_error(SDFCoefficients(returns, factors, include_standard_errors = TRUE, hac_prewhite = TRUE))

  # Testing error handling for incorrect dimensions (transposed matrices).
  expect_error(SDFCoefficients(t(returns), factors, include_standard_errors = TRUE))
  expect_error(SDFCoefficients(returns, t(factors), include_standard_errors = TRUE))
  expect_error(SDFCoefficients(t(returns), t(factors), include_standard_errors = TRUE))

  # Testing errors for wrong input types
  expect_error(SDFCoefficients(c(), factors, include_standard_errors = TRUE))
  expect_error(SDFCoefficients(returns, c(), include_standard_errors = TRUE, hac_prewhite = "c"))
  expect_error(SDFCoefficients(returns, factors, include_standard_errors = "c"))
  expect_error(SDFCoefficients(returns, factors, include_standard_errors = TRUE, hac_prewhite = "c"))

  # Test if the function correctly throws an error when 'returns' has fewer rows than 'factors'.
  expect_error(SDFCoefficients(returns[1:(nrow(returns)-5),], factors))

  # Test if the function correctly throws an error when 'factors' has fewer rows than 'returns'.
  expect_error(SDFCoefficients(returns, factors[1:(nrow(factors)-5),]))

  # Getting results from SDFCoefficients for further validations.
  sdf_coefficients = SDFCoefficients(returns, factors, include_standard_errors = TRUE)

  # Validating the length of the SDF coefficients and standard errors vectors.
  expect_length(sdf_coefficients$sdf_coefficients, n_factors)
  expect_length(sdf_coefficients$standard_errors, n_factors)

  # Comparing computed SDF coefficients with the expected values from manual calculations.
  expect_equal(sdf_coefficients$sdf_coefficients, unname(sdf_coeffs))

})
