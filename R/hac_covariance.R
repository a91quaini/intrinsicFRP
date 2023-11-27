# Author: Alberto Quaini

###############################
######  HAC covariance ########
###############################


#' @title Heteroskedasticity and Autocorrelation robust covariance estimator
#'
#' @name HACcovariance
#' @description This function estimates the long-run covariance matrix of a multivariate
#' centred time series accounting for heteroskedasticity and autocorrelation
#' using the Newey-West (1994)
#' <doi:10.2307/2297912> estimator.
#' The number is selected using the Newey-West plug-in procedure, where
#' `n_lags = 4 * (n_observations/100)^(2/9)`.
#' The function allows to internally prewhiten the series by fitting a VAR(1).
#' All the details can be found in Newey-West (1994)
#' <doi:10.2307/2297912>.
#'
#' @param series A matrix (or vector) of data where each column is a time series.
#' @param prewhite A boolean indicating if the series needs prewhitening by
#' fitting an AR(1). Default is `FALSE`
#' @param check_arguments A boolean `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A symmetric matrix (or a scalar if only one column series is provided)
#' representing the estimated HAC covariance.
#'
#' @examples
#' # Import package data on 6 risk factors and 42 test asset excess returns
#' returns = intrinsicFRP::returns[,-1]
#' factors = intrinsicFRP::factors[,-1]
#'
#' # Fit a linear model of returns on factors
#' fit = stats::lm(returns ~ factors)
#'
#' # Extract residuals from the model
#' residuals = stats::residuals(fit)
#'
#' # Compute the HAC covariance of the residuals
#' hac_covariance = HACcovariance(residuals)
#'
#' # Compute the HAC covariance of the residuals imposing prewhitening
#' hac_covariance_pw = HACcovariance(residuals, prewhite = TRUE)
#'
#' @export
HACcovariance = function(
  series,
  prewhite = FALSE,
  check_arguments = TRUE
) {

  # Check the function arguments if check_arguments is TRUE.
  if (check_arguments) {

    stopifnot("`series` must contain numeric values" = is.numeric(series))
    stopifnot("`series` contains more assets (columns) than observations (rows)" = nrow(series) > ncol(series))
    stopifnot("`series` must not contain missing values (NA/NaN)" = !anyNA(series))
    stopifnot("`prewhite` must be boolean" = is.logical(prewhite))

  }

  # Call the appropriate C++ implementation of the HAC covariance estimation
  # depending whether a vector or a matrix is supplied.
  if (ncol(series) == 1) {

    # HAC variance of a scalar series.
    return(.Call(`_intrinsicFRP_HACVarianceCpp`,
      series,
      prewhite
    ))

  }

  # HAC covariance of multiple series.
  return(.Call(`_intrinsicFRP_HACCovarianceMatrixCpp`,
    series,
    prewhite
  ))

}
