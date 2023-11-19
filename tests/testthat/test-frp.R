test_that("Test FRP", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  # Calculating factor loadings using Fama-MacBeth procedure.
  factor_loadings = t(solve(stats::cov(factors), stats::cov(factors, returns)))
  variance_returns = stats::cov(returns)
  mean_returns = matrix(colMeans(returns), n_returns, 1)

  # Computing risk premia using the standard Fama-MacBeth approach.
  risk_premia = solve(
    t(factor_loadings) %*% factor_loadings,
    t(factor_loadings) %*% mean_returns
  )

  # Calculating risk premia using the Kan-Robotti-Shanken approach.
  var_ret_inv_factor_loadings = solve(variance_returns, factor_loadings)
  krs_risk_premia = solve(
    t(factor_loadings) %*% var_ret_inv_factor_loadings,
    t(var_ret_inv_factor_loadings) %*% mean_returns
  )

  # Testing basic functionality of FRP without errors.
  expect_no_error(FRP(returns, factors))

  # Testing FRP with standard errors included.
  expect_no_error(FRP(returns, factors, include_standard_errors = TRUE))

  # Testing FRP without misspecification robustness but including standard errors.
  expect_no_error(FRP(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE))

  # Testing error handling for incorrect dimensions (transposed matrices).
  expect_error(FRP(t(returns), factors, include_standard_errors = TRUE))
  expect_error(FRP(returns, t(factors), include_standard_errors = TRUE))
  expect_error(FRP(t(returns), t(factors), include_standard_errors = TRUE))

  # Getting results from FRP for further validations.
  krs_frp = FRP(returns, factors, include_standard_errors = TRUE)
  frp = FRP(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE)

  # Validating the length of the risk premia and standard errors vectors.
  expect_length(frp$risk_premia, n_factors)
  expect_length(frp$standard_errors, n_factors)

  expect_length(krs_frp$risk_premia, n_factors)
  expect_length(krs_frp$standard_errors, n_factors)

  # Comparing computed risk premia with the expected values from manual calculations.
  expect_equal(frp$risk_premia, unname(risk_premia))
  expect_equal(krs_frp$risk_premia, unname(krs_risk_premia))

})
