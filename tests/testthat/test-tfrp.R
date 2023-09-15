test_that("Test TFRP", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  covariance_factors_returns = stats::cov(factors, returns)
  variance_returns = stats::cov(returns)
  mean_returns = matrix(colMeans(returns), n_returns, 1)

  risk_premia = covariance_factors_returns %*% solve(
    variance_returns, mean_returns
  )
  risk_premia = unname(risk_premia)

  expect_no_error(TFRP(returns, factors, include_standard_errors = TRUE))
  expect_error(TFRP(t(returns), factors, include_standard_errors = TRUE))
  expect_error(TFRP(returns, t(factors), include_standard_errors = TRUE))
  expect_error(TFRP(t(returns), t(factors), include_standard_errors = TRUE))

  ifrp = TFRP(returns, factors, include_standard_errors = TRUE)

  expect_length(ifrp$risk_premia, n_factors)
  expect_length(ifrp$standard_errors, n_factors)

  expect_equal(ifrp$risk_premia, risk_premia, tolerance = 1e-8)
  expect_equal(
    ifrp$risk_premia,
    TFRP(returns, factors, include_standard_errors = FALSE)$risk_premia,
    tolerance = 1e-8
  )

})
