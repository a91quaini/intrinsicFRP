test_that("Test FRP", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  factor_loadings = t(solve(stats::cov(factors), stats::cov(factors, returns)))
  variance_returns = stats::cov(returns)
  mean_returns = matrix(colMeans(returns), n_returns, 1)

  risk_premia = solve(
    t(factor_loadings) %*% factor_loadings,
    t(factor_loadings) %*% mean_returns
  )

  var_ret_inv_factor_loadings = solve(variance_returns, factor_loadings)
  krs_risk_premia = solve(
    t(factor_loadings) %*% var_ret_inv_factor_loadings,
    t(var_ret_inv_factor_loadings) %*% mean_returns
  )

  expect_no_error(FRP(returns, factors, include_standard_errors = TRUE))
  expect_no_error(FRP(returns, factors, krs = FALSE, include_standard_errors = TRUE))

  expect_error(FRP(t(returns), factors, include_standard_errors = TRUE))
  expect_error(FRP(returns, t(factors), include_standard_errors = TRUE))
  expect_error(FRP(t(returns), t(factors), include_standard_errors = TRUE))

  krs_frp = FRP(returns, factors, include_standard_errors = TRUE)
  frp = FRP(returns, factors, krs = FALSE, include_standard_errors = TRUE)

  expect_length(frp$risk_premia, n_factors)
  expect_length(frp$standard_errors, n_factors)

  expect_length(krs_frp$risk_premia, n_factors)
  expect_length(krs_frp$standard_errors, n_factors)

  expect_equal(frp$risk_premia, unname(risk_premia))
  expect_equal(krs_frp$risk_premia, unname(krs_risk_premia))

})
