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
// using the Newey-West estimator.
// If the number of lags is not provided, they are selected using the Newey-West (1994)
// <doi:10.2307/2297912> plug-in procedure, where
// `n_lags = 4 * (n_observations/100)^(2/9)`.
// The function allows to internally prewhiten the series by fitting a VAR(1).
//
// @param series A matrix of data where each column is a time series of model residuals.
// @param n_lags An integer indicating the number of lags. Default is -1.
// @param prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1). Default is false.
//
// @return A symmetric matrix representing the estimated HAC covariance matrix.
// [[Rcpp::export]]
arma::mat HACCovarianceMatrixCpp(
  arma::mat& series,
  int n_lags = -1,
  const bool prewhite = false
);

// Function for internal use
//
// Computes the standard errors from the HAC covariance matrix.
//
// This function derives the HAC standard errors for each series by taking the square root
// of the diagonal elements of the HAC covariance matrix. These standard errors are crucial
// for statistical inference when addressing heteroskedasticity and autocorrelation.
// If the number of lags is not provided, they are selected using the Newey-West (1994) <doi:10.2307/2297912>
// plug-in procedure, where `n_lags = 4 * (n_observations/100)^(2/9)`.
//
// @param series A matrix of data where each column is a time series of model residuals.
// @param n_lags An integer indicating the number of lags. Default is -1.
//
// @return A vector of HAC marginal standard errors for each series.
arma::vec HACStandardErrorsCpp(
  const arma::mat& series,
  int n_lags = -1
);

// Function for internal use
//
// Computes the Heteroskedasticity and Autocorrelation robust variance.
//
// This function calculates the variance of a univariate time series,
// adjusting for autocorrelation
// and potential heteroskedasticity using the Newey-West estimator.
// If the number of lags is not provided, they are selected using the Newey-West (1994) <doi:10.2307/2297912>
// plug-in procedure, where `n_lags = 4 * (n_observations/100)^(2/9)`.
// The function allows to internally prewhiten the series by fitting a VAR(1).
//
// @param series A matrix of data where each column is a time series of model residuals.
// @param n_lags An integer indicating the number of lags. Default is -1.
// @param prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1). Default is false.
//
// @return The Newey-West adjusted variance of the input series.
// [[Rcpp::export]]
double HACVarianceCpp(
  arma::vec& series,
  int n_lags = -1,
  const bool prewhite = false
);

///////////////////////
//// Prewhitening /////

// Function for internal use
// pre-whitening by fitting a VAR(1)
void HACPrewhiteningCpp(arma::mat& series, arma::mat& coefficients);

// Function for internal use
// pre-whitening by fitting a VAR(1)
void HACPrewhiteningScalarCpp(arma::vec& series, double coefficient);

// Function for internal use
// obtain HAC covariance estimation of pre-whitened series
void HACRevertPrewhiteningCpp(
  const arma::mat& coefficients,
  arma::mat& hac_covariance
);

// Function for internal use
// obtain HAC variance estimation of pre-whitened series
void HACRevertPrewhiteningScalarCpp(
  const double coefficient,
  double hac_covariance
);

#endif
