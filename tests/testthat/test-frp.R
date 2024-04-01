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

  # Test if prewhite works
  expect_no_error(FRP(returns, factors, include_standard_errors = TRUE, hac_prewhite = TRUE))

  # Testing FRP without misspecification robustness but including standard errors.
  expect_no_error(FRP(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE))
  expect_no_error(FRP(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE, target_level_gkr2014_screening = 0.05))
  expect_no_error(FRP(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE, target_level_gkr2014_screening = 0.95))


  # Testing GKR factor screening with misspecification_robust
  expect_no_error(FRP(returns, factors, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))
  expect_no_error(FRP(returns, factors, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.95))
  expect_no_error(FRP(returns, factors, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))

  # Testing GKR factor screening with misspecification_robust on simulated factors
  factors_usl = matrix(stats::rnorm(nrow(returns) * 2), nrow(returns), 2)
  expect_no_error(FRP(returns, factors_usl, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))
  expect_no_error(FRP(returns, factors_usl, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.95))
  expect_no_error(FRP(returns, factors_usl, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))


  # Testing error handling for incorrect dimensions (transposed matrices).
  expect_error(FRP(t(returns), factors, include_standard_errors = TRUE))
  expect_error(FRP(returns, t(factors), include_standard_errors = TRUE))
  expect_error(FRP(t(returns), t(factors), include_standard_errors = TRUE))

  # Testing errors for wrong input types
  expect_error(FRP(c(), factors, include_standard_errors = TRUE))
  expect_error(FRP(returns, c(), include_standard_errors = TRUE, hac_prewhite = "c"))
  expect_error(FRP(returns, factors, include_standard_errors = "c"))
  expect_error(FRP(returns, factors, include_standard_errors = TRUE, hac_prewhite = "c"))

  # Test if the function correctly throws an error when 'returns' has fewer rows than 'factors'.
  expect_error(FRP(returns[1:(nrow(returns)-5),], factors))

  # Test if the function correctly throws an error when 'factors' has fewer rows than 'returns'.
  expect_error(FRP(returns, factors[1:(nrow(factors)-5),]))

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
