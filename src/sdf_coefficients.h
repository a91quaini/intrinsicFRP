// Author: Alberto Quaini

#ifndef SDF_COEFFICIENTS_H
#define SDF_COEFFICIENTS_H

#include <RcppArmadillo.h>

///////////////////////////////////////////////
////////////  SDF Coefficients ////////////////
///////////////////////////////////////////////

// Compute SDF Coefficients
//
// @name SDFCoefficientsCpp
// @description Computes the misspecification-robust SDF coefficients of
// Gospodinov-Kan-Robotti (2014) <https://doi.org/10.1093/rfs/hht135>:
// `GKRSDFcoefficients = (C' * V[R]^{-1} * C)^{-1} * C' * V[R]^{-1} * E[R]`
// from data on factors `F` and test
// asset excess returns `R`.
// These notions of SDF coefficients minimize pricing errors:
// `argmin_{d} (E[R] - Cov[R,F] * d)' * V[R]^{-1} * (E[R] - Cov[R,F] * d)`.
// Optionally computes the corresponding
// heteroskedasticity and autocorrelation robust standard errors using the
// Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the
// number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
//
// @param returns A `n_observations x n_returns`-dimensional matrix of test asset
// excess returns.
// @param factors A `n_observations x n_factors`-dimensional matrix of factors.
// @param misspecification_robust A boolean: `TRUE` for the
// "misspecification-robust" Kan-Robotti-Shanken (2013) GLS approach using the
// inverse covariance matrix of returns; `FALSE` for standard Fama-MacBeth
// risk premia. Default is `TRUE`.
// @param include_standard_errors A boolean: `TRUE` if you want to compute the
// SDF coefficients' HAC standard errors; `FALSE` otherwise.
// Default is `FALSE`.
// @param hac_prewhite A boolean indicating if the series needs pre-whitening by
// fitting an AR(1) in the internal heteroskedasticity and autocorrelation
// robust covariance (HAC) estimation. Default is `false`.
// @param target_level_gkr2014_screening A number indicating the target level of
// the tests underlying the factor screening procedure in Gospodinov-Kan-Robotti
// (2014). If it is zero, then no factor screening procedure is
// implemented. Otherwise, it implements an iterative screening procedure
// based on the sequential removal of factors associated with the smallest insignificant
// t-test of a nonzero SDF coefficient. The threshold for the absolute t-test is
// `target_level_gkr2014_screening / n_factors`, where n_factors indicate the
// number of factors in the model at the current iteration. Default is `0.`, i.e.,
// no factor screening.
//
// @return A list containing `n_factors`-dimensional vector of fSDF coefficients
// in `"sdf_coefficients"`; if `include_standard_errors = TRUE`, then
// it further includes `n_factors`-dimensional vector of SDF coefficients'
// standard errors in `"standard_errors"`;
// if `target_level_gkr2014_screening >= 0`, it further includes the indices of
// the selected factors in `selected_factor_indices`.
//
// [[Rcpp::export]]
Rcpp::List SDFCoefficientsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust = true,
  const bool include_standard_errors = false,
  const bool hac_prewhite = false,
  const double target_level_gkr2014_screening = 0.
);

// Function for internal use
// Compiles the list returned by `FRPCpp`.
Rcpp::List ReturnSDFCoefficientsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust,
  const bool include_standard_errors,
  const bool hac_prewhite = false
);

// Function for internal use
// Computes the Fama-MacBeth (1973) <doi:10.1086/260061> from moments
// extracted from data
arma::vec FMSDFCoefficientsCpp(
  const arma::mat& covariance_returns_factors,
  const arma::vec& mean_returns
);

// Function for internal use
// Computes the misspecification robust SDF coefficients of
// Gospodinov-Kan-Robotti (2014)
// <https://doi.org/10.1093/rfs/hht135> from moments extracted from data
arma::vec GKRSDFCoefficientsCpp(
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
);

// Function for internal use
// Compute the standard errors of the Fama-MacBeth (1973) <doi:10.1086/260061>
// SDF coefficients from moments using the heteroskedasticity and
// autocorrelation robust standard errors using the Newey-West (1994)
// <doi:10.2307/2297912> plug-in procedure to select the number of relevant lags,
// i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
// Formulas adapted for working with excess returns are own derived.
arma::vec StandardErrorsFMSDFCoefficientsCpp(
  const arma::vec& fm_sdf_coefficients,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::vec& mean_returns,
  const bool hac_prewhite = false
);

// Function for internal use
// Compute the standard errors of the Gospodinov-Kan-Robotti (2014) <https://doi.org/10.1093/rfs/hht135>
// SDF coefficients from moments using the heteroskedasticity and
// autocorrelation robust standard errors using the Newey-West (1994)
// <doi:10.2307/2297912> plug-in procedure to select the number of relevant lags,
// i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
// Formula adapted for using excess returns are given in Kan Robotti (2008)
// <https://doi.org/10.1016/j.jempfin.2008.03.003>
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
