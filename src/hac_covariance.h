// Author: Alberto Quaini

#ifndef HAC_COVARIANCE_H
#define HAC_COVARIANCE_H

#include <RcppArmadillo.h>

// Function for internal use
//
// Computes the heteroskedasticity and autocorrelation consistent (HAC) covariance matrix.
//
// This function estimates the long-run covariance matrix of a multivariate
// centred time series accounting for heteroskedasticity and autocorrelation
// using the Newey-West (1994)
// <doi:10.2307/2297912> estimator.
// The number is selected using the Newey-West plug-in procedure, where
// `n_lags = 4 * (n_observations/100)^(2/9)`.
// The function allows to internally prewhiten the series by fitting a VAR(1).
// All the details can be found in Newey-West (1994)
// <doi:10.2307/2297912>.
//
// @param series A matrix (or vector) of data where each column is a time series.
// @param prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1). Default is `false`.
//
// @return A symmetric matrix representing the estimated HAC covariance matrix.
// [[Rcpp::export]]
arma::mat HACCovarianceMatrixCpp(
  arma::mat& series,
  const bool prewhite = false
);

// Function for internal use
//
// Computes the standard errors from the HAC covariance matrix.
//
// This function derives the HAC standard errors for each series by taking the square root
// of the diagonal elements of the HAC covariance matrix. These standard errors are crucial
// for statistical inference when addressing heteroskedasticity and autocorrelation.
// The number of lags is selected using the Newey-West (1994) <doi:10.2307/2297912>
// plug-in procedure, where `n_lags = 4 * (n_observations/100)^(2/9)`.
//
// @param series A matrix of data where each column is a time series of model residuals.
// @param prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1). Default is `false`.
//
// @return A vector of HAC marginal standard errors for each series.
arma::vec HACStandardErrorsCpp(
  arma::mat& series,
  const bool prewhite = false
);

// Function for internal use
//
// Computes the Heteroskedasticity and Autocorrelation robust variance
// of a scalar series.
//
// This function calculates the variance of a univariate time series,
// adjusting for autocorrelation
// and potential heteroskedasticity using the Newey-West estimator.
// The number is selected using the Newey-West (1994) <doi:10.2307/2297912>
// plug-in procedure, where `n_lags = 4 * (n_observations/100)^(2/9)`.
// The function allows to internally prewhiten the series by fitting a VAR(1).
//
// @param series A matrix of data where each column is a time series of model residuals.
// @param n_lags An integer indicating the number of lags. If it is negative,
// then `n_lags = 4 * (n_observations/100)^(2/9)`. Default is `-1`.
// @param prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1). Default is `false`.
//
// @return The Newey-West adjusted variance of the input series.
// [[Rcpp::export]]
double HACVarianceCpp(
  arma::vec& series,
  const bool prewhite = false
);

///////////////////////
//// Prewhitening /////

// Function for internal use
// pre-whitening of matrix series by fitting a VAR(1)
void HACPrewhiteningCpp(arma::mat& series, arma::mat& coefficients);

// Function for internal use
// marginal pre-whitening of matrix series by fitting a VAR(1)
void HACPrewhiteningCpp(arma::mat& series, arma::vec& coefficients);

// Function for internal use
// pre-whitening of scalar series by fitting a VAR(1)
void HACPrewhiteningCpp(arma::vec& series, double coefficient);

// Function for internal use
// obtain HAC covariance estimation of pre-whitened matrix series
void HACRevertPrewhiteningCpp(
  const arma::mat& coefficients,
  arma::mat& hac_covariance
);

// Function for internal use
// obtain the diagonal HAC covariance estimation of pre-whitened matrix series
void HACRevertPrewhiteningCpp(
  const arma::vec& coefficients,
  arma::rowvec& hac_covariance
);

// Function for internal use
// obtain HAC variance estimation of pre-whitened scalar series
void HACRevertPrewhiteningCpp(
  const double coefficient,
  double hac_covariance
);

#endif
