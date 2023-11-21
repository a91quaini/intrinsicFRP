test_that("Test HACcovariance", {

  returns <- intrinsicFRP::returns[,-1]
  factors <- intrinsicFRP::factors[,-1]
  fit <- stats::lm(returns ~ factors)
  residuals <- stats::residuals(fit)

  hac_covariance <- HACcovariance(residuals)

  # HACcovariance returns correct matrix dimensions
  expect_true(is.matrix(hac_covariance))
  expect_equal(ncol(hac_covariance), ncol(residuals))
  expect_equal(nrow(hac_covariance), ncol(residuals))

  # HACcovariance handles single column series
  single_series <- residuals[,1, drop = FALSE]

  hac_covariance_single <- HACcovariance(single_series)
  expect_true(is.numeric(hac_covariance_single))
  expect_length(hac_covariance_single, 1)

  # HACcovariance handles prewhitening
  hac_covariance_pw <- HACcovariance(residuals, prewhite = TRUE)
  expect_true(is.matrix(hac_covariance_pw))
  expect_equal(dim(hac_covariance_pw), c(ncol(residuals), ncol(residuals)))

  hac_covariance_pw_single <- HACcovariance(single_series, prewhite = TRUE)
  expect_true(is.numeric(hac_covariance_pw_single))
  expect_length(hac_covariance_pw_single, 1)

  # HACcovariance errors with invalid inputs
  expect_error(HACcovariance(NULL))
  expect_error(HACcovariance(matrix(NA, nrow = 10, ncol = 5)))
  expect_error(HACcovariance(residuals, prewhite = "yes"))

  # HACcovariance for n_observations smaller than 5
  expect_error(HACcovariance(residuals[1:3,]))
  expect_error(HACcovariance(single_series[1:3, drop = FALSE]))

})
