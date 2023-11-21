// Author: Alberto Quaini

#ifndef GKR_FACTOR_SCREENING_H
#define GKR_FACTOR_SCREENING_H

#include <RcppArmadillo.h>

// Perform the factor screening procedure of Gospodinov-Kan-Robotti (2014)
// from moments extracted from data
//
// @name GKRFactorScreeningCpp
// @description Performs the factor screening procedure of
// Gospodinov-Kan-Robotti (2014) <doi:10.2139/ssrn.2579821>, which is
// an iterative screening procedure
// based on the sequential removal of factors associated with the smallest insignificant
// t-test of a nonzero misspecification-robust SDF coefficient. The significance threshold for the
// absolute t-test is given by `target_level_gkr2014_screening / n_factors`,
// where n_factors indicate the number of factors in the model at the current iteration;
// that is, it takes care of the multiple testing problem via a conservative
// Bonferroni correction. Standard errors are computed with the
// heteroskedasticity and autocorrelation using the Newey-West estimator.
// The number is selected using the Newey-West (1994)
// <doi:10.2307/2297912> plug-in procedure, where
// `n_lags = 4 * (n_observations/100)^(2/9)`.
// The function allows to internally prewhiten the series by fitting a VAR(1).
//
// @param returns `n_observations x n_returns`-dimensional matrix of test asset
// excess returns.
// @param factors `n_observations x n_factors`-dimensional matrix of risk
// factors.
// @param target_level Double specifying the target level used for rank estimation.
// @param hac_prewhite A boolean indicating if the series needs prewhitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// robust covariance (HAC) estimation. Default is `false`.
//
// @return A list contaning the GKR SDF coefficients, their standard errors and
// squared t-statistics. It further contain the indices of the selected factors.
//
// [[Rcpp::export]]
Rcpp::List GKRFactorScreeningCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level = 0.05,
  const bool hac_prewhite = false
);

// Function for internal use
// It Performs the factor screening procedure of
// Gospodinov-Kan-Robotti (2014) <doi:10.2139/ssrn.2579821>
// based on moments computed on data comprising risk factors and test asset
// excess returns.
Rcpp::List GKRFactorScreeningCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double target_level = 0.05,
  const bool hac_prewhite = false
);

// Function for internal use
// Computes the misspecification-robust SDF coefficients of
// Gospodinov-Kan-Robotti (2014) <doi:10.2139/ssrn.2579821>
// from moments computed on data comprising risk factors and test asset
// excess returns.
arma::vec GKRSDFCoefficientsCpp(
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// Function for internal use
// Compute the standard errors of the Gospodinov-Kan-Robotti (2014) <doi:10.2139/ssrn.2579821>
// misspecification-robust SDF coefficients !multiplied by sqrt(n_observations)!
// from moments computed on data comprising
// risk factors and test asset excess returns. It uses the heteroskedasticity and
// autocorrelation robust standard errors using the Newey-West (1994)
// <doi:10.2307/2297912> plug-in procedure to select the number of relevant lags,
// i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
// If the series needs prewhitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// The function allows to internally prewhiten the series by fitting a VAR(1).
arma::vec StandardErrorsGKRSDFCoefficientsCpp(
  const arma::vec& gkr_sdf_coefficients,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite = false
);

#endif
