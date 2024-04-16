# Author: Alberto Quaini

##############################################################
####  Test of new factors by Feng Giglio and Xiu (2020) ######
##############################################################

#' @title Testing for the pricing contribution of new factors.
#'
#' @name FGXFactorsTest
#' @description Computes the three-pass procedure of Feng Giglio and Xiu (2020)
#' <doi:org/10.1111/jofi.12883>, which evaluates the contribution
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
#' @param gross_returns A `n_observations x n_returns`-dimensional matrix of test asset
#' gross returns.
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
#' `"standard_errors"`; it also returns the index of control factors that are
#' selected by the two-step selection procedure.
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
#' control_factors = intrinsicFRP::factors[,2:4]
#' new_factors = intrinsicFRP::factors[,5:7]
#' returns = intrinsicFRP::returns[,-1]
#' RF = intrinsicFRP::risk_free[,-1]
#'
#' gross_returns = returns + 1 + RF
#'
#' output = FGXFactorsTest(
#'   gross_returns,
#'   control_factors,
#'   new_factors
#' )
#'
#' @export
FGXFactorsTest = function(
  gross_returns,
  control_factors,
  new_factors,
  n_folds = 5,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {

    CheckData(gross_returns, control_factors)
    stopifnot("`new_factors` must contain numeric values" = is.numeric(new_factors))
    stopifnot("`new_factors` contains more variables (columns) than observations (rows)" = nrow(new_factors) > ncol(new_factors))
    stopifnot("`new_factors` must not contain missing values (NA/NaN)" = !anyNA(new_factors))
    stopifnot("`n_folds` must be numeric" = is.numeric(n_folds))

  }

  # set number of new factors
  n_new = ncol(new_factors)

  # compute moments of data
  avg_returns = colMeans(gross_returns)
  cov_returns_control_factors = stats::cov(gross_returns, control_factors)
  cov_returns_new_factors = stats::cov(gross_returns, new_factors)

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
  for (ii in 1:n_new) {

    second_lasso = glmnet::cv.glmnet(
      x = cov_returns_control_factors,
      y = cov_returns_new_factors[, ii],
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
  predictors = cbind(
    matrix(1, ncol(gross_returns), 1),
    cov_returns_new_factors,
    cov_returns_control_factors[, idx_selected, drop=FALSE]
  )

  # cross_sectional_fit = stats::lm(avg_returns ~ covariances)
  # compute the SDF coefficients
  sdf_coefficients = solve(
    t(predictors) %*% predictors,
    t(predictors) %*% avg_returns
  )

  ## Estimate the Newey-West type covariance estimator of the new factors'
  ## SDF coefficients

  # if no control factor was selected
  if (length(idx_selected) == 0) {

    # compute the Newey-West type covariance estimator
    covariance = .Call(`_intrinsicFRP_FGXThreePassCovarianceNoControlsCpp`,
      gross_returns,
      new_factors,
      sdf_coefficients[-1]
    )

  # if once control factor was selected
  } else if (length(idx_selected) == 1) {

    # compute the Newey-West type covariance estimator
    covariance = .Call(`_intrinsicFRP_FGXThreePassCovarianceCpp`,
      gross_returns,
      control_factors[, idx_selected, drop = FALSE],
      new_factors,
      sdf_coefficients[-1]
    )

  } else {

    # lasso selection:
    # lasso regression of each new factors on the control factors.
    control_factors = control_factors[,idx_selected]
    cov_selected = c()

    for (ii in 1:n_new) {

      lasso_fit = glmnet::cv.glmnet(
        x = control_factors,
        y = new_factors[, ii],
        n_folds = n_folds
      )

      lasso_fit_coeffs = stats::coef(lasso_fit, s = "lambda.min")[-1]
      cov_selected = union(
        cov_selected,
        which(lasso_fit_coeffs != 0, arr.ind = T)
      )

    }

    if (length(cov_selected) == 0) {

      # compute the Newey-West type covariance estimator
      covariance = .Call(`_intrinsicFRP_FGXThreePassCovarianceNoControlsCpp`,
        gross_returns,
        new_factors,
        sdf_coefficients[1 + 1:n_new]
      )

    } else {

      # index selected sdf_coefficients
      idx_sdf_coeff = 1 + c(1:n_new, cov_selected + n_new)

      # compute the Newey-West type covariance estimator
      covariance = .Call(`_intrinsicFRP_FGXThreePassCovarianceCpp`,
        gross_returns,
        control_factors[, cov_selected, drop = FALSE],
        new_factors,
        sdf_coefficients[idx_sdf_coeff]
      )

    }

  }

  # extract the standard errors
  standard_errors = sqrt(diag(covariance)) / sqrt(nrow(gross_returns))

  # return a list containing the new factors' SDF coefficients and corresponding
  # standard errors
  output = list(sdf_coefficients[1 + 1:n_new], standard_errors, idx_selected)
  names(output) = c("sdf_coefficients", "standard_errors", "controls_selected")

  return(output)

}
