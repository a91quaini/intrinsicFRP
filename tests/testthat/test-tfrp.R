test_that("Test TFRP", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  # Calculating the covariance between factors and returns, and the variance of returns.
  covariance_factors_returns = stats::cov(factors, returns)
  variance_returns = stats::cov(returns)
  mean_returns = matrix(colMeans(returns), n_returns, 1)

  # Computing risk premia using the tradable factor risk premia formula.
  risk_premia = covariance_factors_returns %*% solve(variance_returns, mean_returns)
  risk_premia = unname(risk_premia)

  # Testing basic functionality of TFRP without errors including standard errors.
  expect_no_error(TFRP(returns, factors, include_standard_errors = TRUE))

  # Testing error handling for incorrect dimensions (transposed matrices).
  expect_error(TFRP(t(returns), factors, include_standard_errors = TRUE))
  expect_error(TFRP(returns, t(factors), include_standard_errors = TRUE))
  expect_error(TFRP(t(returns), t(factors), include_standard_errors = TRUE))

  # Getting results from TFRP for further validations.
  ifrp = TFRP(returns, factors, include_standard_errors = TRUE)

  # Validating the length of the risk premia and standard errors vectors.
  expect_length(ifrp$risk_premia, n_factors)
  expect_length(ifrp$standard_errors, n_factors)

  # Comparing computed risk premia with the expected values from manual calculations.
  expect_equal(ifrp$risk_premia, risk_premia, tolerance = 1e-8)

  # Ensuring consistency in the risk premia calculation with and without standard errors.
  expect_equal(
    ifrp$risk_premia,
    TFRP(returns, factors, include_standard_errors = FALSE)$risk_premia,
    tolerance = 1e-8
  )

})
