# Author: Alberto Quaini

##############################################################
####  Test of new factors by Feng Giglio and Xiu (2020) ######
##############################################################

#' @title Testing for the pricing contribution of new factors.
#'
#' @name FGXFactorsTest
#' @description Computes the three-pass procedure of Feng Giglio and Xiu (2020)
#' <doi:https://doi.org/10.1111/jofi.12883>, which evaluates the contribution
#' to cross-sectional pricing of any new factors on top of a set of control
#' factors.
#' The third step is a OLS regression of average returns on the covariances between
#' asset returns and the new factors, as well as the control factors selected
#' in either one of the first two steps.
#' The first stwo steps consists in (i) a Lasso regression of average returns
#' on the ovariances between asset returns and all control factors and (ii)
#' a Lasso regression of the covariances between asset returns and the new factors
#' on the ovariances between asset returns and all control factors.
#' The second selection aims at correcting for potential omitted variables in the
#' first selection.
#' Tuning of the penalty parameters in the Lasso regressions is performed via
#' Cross Validation (CV).
#' Standard errors are computed following Feng Giglio and Xiu (2020) using the
#' Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the number
#' of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#' For the standard error computations, the function allows to internally
#' pre-whiten the series by fitting a VAR(1),
#' i.e., a vector autoregressive model of order 1.
#'
#' @param returns A `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param control_factors A `n_observations x n_control_factors`-dimensional
#' matrix of control or benchmark factors.
#' @param new_factors A `n_observations x n_new_factors`-dimensional
#' matrix of new factors.
#' @param n_folds An integer indicating
#' the number of k-fold for cross validation. Default is `5`.
#' @param check_arguments A boolean `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A list containing the `n_new_factors`-dimensional vector of SDF
#' coefficients in `"sdf_coefficients"` and corresponding standard errors in
#' `"standard_errors"`.
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
#' control_factors = intrinsicFRP::factors[,2:4]
#' new_factors = intrinsicFRP::factors[,5:7]
#' returns = intrinsicFRP::returns[,-1]
#'
#' output = FGXFactorsTest(
#'   returns,
#'   control_factors,
#'   new_factors
#' )
#'
#' @export
FGXFactorsTest = function(
  returns,
  control_factors,
  new_factors,
  n_folds = 5,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {

    CheckData(returns, control_factors)
    stopifnot("`new_factors` must contain numeric values" = is.numeric(new_factors))
    stopifnot("`new_factors` contains more variables (columns) than observations (rows)" = nrow(new_factors) > ncol(new_factors))
    stopifnot("`new_factors` must not contain missing values (NA/NaN)" = !anyNA(new_factors))
    stopifnot("`n_folds` must be numeric" = is.numeric(n_folds))

  }

  # compute moments of data
  avg_returns = colMeans(returns)
  cov_returns_control_factors = stats::cov(returns, control_factors)
  cov_returns_new_factors = stats::cov(returns, new_factors)

  # first lasso selection:
  # lasso regression of average returns on the covariances between returns
  # and all control factors
  first_lasso = glmnet::cv.glmnet(
    x = cov_returns_control_factors,
    y = avg_returns,
    n_folds = n_folds
  )

  first_lasso_coeffs = stats::coef(first_lasso, s = "lambda.min")[-1]
  idx_selected = which(first_lasso_coeffs != 0, arr.ind = T)

  # second lasso selection:
  # lasso regression of each covariance between new factors and returns on
  # covariances between control factors and returns
  # take the union of the first and all second step selections
  for (ii in 1:ncol(new_factors)) {

    second_lasso = glmnet::cv.glmnet(
      x = cov_returns_control_factors,
      y = cov_returns_new_factors[, ii, drop=FALSE],
      n_folds = n_folds
    )

    second_lasso_coeffs = stats::coef(second_lasso, s = "lambda.min")[-1]
    idx_selected = union(
      idx_selected,
      which(second_lasso_coeffs != 0, arr.ind = T)
    )

  }

  # second-step cross-sectional regression:
  # OLS regression of average returns on covariances between new factors and returns
  # and on the covariances between the selected control factors and returns.
  # the selected control factors are the union of the first and all second step selections.
  covariances = cbind(
    cov_returns_new_factors,
    cov_returns_control_factors[, idx_selected, drop=FALSE]
  )

  # cross_sectional_fit = stats::lm(avg_returns ~ covariances)
  # compute the SDF coefficients
  sdf_coefficients = solve(
    t(covariances) %*% covariances,
    t(covariances) %*% avg_returns
  )

  ## Estimate the Newey-West type covariance estimator of the new factors'
  ## SDF coefficients

  # lasso selection:
  # lasso regression of each new factors on the control factors.
  control_factors = control_factors[,idx_selected]
  cov_selected = c()

  for (ii in 1:ncol(new_factors)) {

    lasso_fit = glmnet::cv.glmnet(
      x = control_factors,
      y = new_factors[, ii, drop=FALSE],
      n_folds = n_folds
    )

    lasso_fit_coeffs = stats::coef(lasso_fit, s = "lambda.min")[-1]
    cov_selected = union(
      cov_selected,
      which(lasso_fit_coeffs != 0, arr.ind = T)
    )

  }

  # turn R indexing into cpp indexing
  cov_selected = cov_selected - 1

  # compute the Newey-West type covariance estimator
  covariance = .Call(`_intrinsicFRP_FGXThreePassCovarianceCpp`,
    returns,
    control_factors,
    new_factors,
    sdf_coefficients,
    cov_selected
  )

  # extract the standard errors
  standard_errors = sqrt(diag(covariance)) / sqrt(nrow(returns))

  # return a list containing the new factors' SDF coefficients and corresponding
  # standard errors
  output = list(sdf_coefficients[1:ncol(new_factors)], standard_errors)
  names(output) = c("sdf_coefficients", "standard_errors")

  return(output)

}
