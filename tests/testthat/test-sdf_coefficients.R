test_that("Test SDF coefficients", {

  factors = factors[,-1]
  returns = returns[,-1]

  n_factors = ncol(factors)
  n_returns = ncol(returns)

  # Calculating factor loadings using Fama-MacBeth procedure.
  covariance_returns_factors = stats::cov(returns, factors)
  variance_returns = stats::cov(returns)
  mean_returns = matrix(colMeans(returns), n_returns, 1)

  # Computing SDF coefficients using the standard Fama-MacBeth approach.
  fm_coeff = solve(
    t(covariance_returns_factors) %*% covariance_returns_factors,
    t(covariance_returns_factors) %*% mean_returns
  )

  # Calculating SDF coefficients using the Gospodinov-Kan-Robotti approach.
  var_ret_inv_cov_ret_fac = solve(variance_returns, covariance_returns_factors)
  gkr_coeff = solve(
    t(var_ret_inv_cov_ret_fac) %*% covariance_returns_factors,
    t(var_ret_inv_cov_ret_fac) %*% mean_returns
  )

  # Testing basic functionality of FRP without errors.
  expect_no_error(SDFCoefficients(returns, factors))

  # Testing SDFCoefficients with standard errors included.
  expect_no_error(SDFCoefficients(returns, factors, include_standard_errors = TRUE))

  # Test if prewhite works
  expect_no_error(SDFCoefficients(returns, factors, include_standard_errors = TRUE, hac_prewhite = TRUE))

  # Testing SDFCoefficients without misspecification robustness but including standard errors.
  expect_no_error(SDFCoefficients(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE))
  expect_no_error(SDFCoefficients(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE, target_level_gkr2014_screening = 0.05))
  expect_no_error(SDFCoefficients(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE, target_level_gkr2014_screening = 0.95))

  # Testing GKR factor screening with misspecification_robust
  expect_no_error(SDFCoefficients(returns, factors, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))
  expect_no_error(SDFCoefficients(returns, factors, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.95))
  expect_no_error(SDFCoefficients(returns, factors, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))

  # Testing GKR factor screening with misspecification_robust for simulated data
  factors_usl = matrix(stats::rnorm(nrow(returns) * 2), nrow(returns), 2)
  expect_no_error(SDFCoefficients(returns, factors_usl, misspecification_robust = FALSE, target_level_gkr2014_screening = 0.95))
  expect_no_error(SDFCoefficients(returns, factors_usl, misspecification_robust = FALSE, target_level_gkr2014_screening = 0.05))
  expect_no_error(SDFCoefficients(returns, factors_usl, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.95))
  expect_no_error(SDFCoefficients(returns, factors_usl, misspecification_robust = TRUE, target_level_gkr2014_screening = 0.05))


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

  # Getting results from FRP for further validations.
  gkr_sdf_coefficients = SDFCoefficients(returns, factors, include_standard_errors = TRUE)
  fm_sdf_coefficients = SDFCoefficients(returns, factors, misspecification_robust = FALSE, include_standard_errors = TRUE)

  # Validating the length of the risk premia and standard errors vectors.
  expect_length(fm_sdf_coefficients$sdf_coefficients, n_factors)
  expect_length(fm_sdf_coefficients$standard_errors, n_factors)

  expect_length(gkr_sdf_coefficients$sdf_coefficients, n_factors)
  expect_length(gkr_sdf_coefficients$standard_errors, n_factors)

  # Comparing computed risk premia with the expected values from manual calculations.
  expect_equal(fm_sdf_coefficients$sdf_coefficients, unname(fm_coeff))
  expect_equal(gkr_sdf_coefficients$sdf_coefficients, unname(gkr_coeff))

})
