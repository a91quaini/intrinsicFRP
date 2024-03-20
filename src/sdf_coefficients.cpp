// Author: Alberto Quaini

#include "sdf_coefficients.h"
#include "hac_covariance.h"
#include "utils.h"

/////////////////////////////
///// SDFCoefficientsCpp ////

Rcpp::List SDFCoefficientsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  // Check if standard errors are to be included in the calculation
  if (include_standard_errors) {

    // compute moments and SDF coefficients
    const arma::mat covariance_returns_factors = arma::cov(returns, factors);
    const arma::mat variance_returns = arma::cov(returns);
    const arma::vec mean_returns = arma::mean(returns).t();
    const arma::vec sdf_coefficients = GKRSDFCoefficientsCpp(
      covariance_returns_factors,
      variance_returns,
      mean_returns
    );

    // return SDF coefficients and their standard errors
    return Rcpp::List::create(
      Rcpp::Named("sdf_coefficients") = sdf_coefficients,
      Rcpp::Named("standard_errors") = StandardErrorsGKRSDFCoefficientsCpp(
        sdf_coefficients,
        returns,
        factors,
        covariance_returns_factors,
        variance_returns,
        mean_returns,
        hac_prewhite
      )
    );

  }

  // return SDF coefficients
  return Rcpp::List::create(
    Rcpp::Named("sdf_coefficients") = GKRSDFCoefficientsCpp(
      arma::cov(returns, factors),
      arma::cov(returns),
      arma::mean(returns).t()
    )
  );

}

////////////////////////////////
///// GKRSDFCoefficientsCpp ////

arma::vec GKRSDFCoefficientsCpp(
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  // Solve for `V[R]^{-1} * Cov[R, F]
  const arma::mat var_ret_inv_cov_ret_fac = SolveSympd(
    variance_returns,
    covariance_returns_factors
  );

  // Compute sdf coefficients using the GKR method
  return SolveSympd(
    var_ret_inv_cov_ret_fac.t() * covariance_returns_factors,
    var_ret_inv_cov_ret_fac.t()
  ) * mean_returns;

}

//////////////////////////////////////////////
///// StandardErrorsGKRSDFCoefficientsCpp ////

arma::vec StandardErrorsGKRSDFCoefficientsCpp(
  const arma::vec& gkr_sdf_coefficients,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Compute the inverse of variance_returns and useful matrices that use it
  const arma::mat var_ret_inv = InvSympd(variance_returns);
  const arma::mat var_ret_inv_cov_ret_fac = var_ret_inv * covariance_returns_factors;
  const arma::vec var_ret_inv_mean_ret = var_ret_inv * mean_returns;

  // Compute matrix H used in the standard error calculation
  const arma::mat h_matrix = InvSympd(
    covariance_returns_factors.t() * var_ret_inv_cov_ret_fac
  );

  // Compute matrix A used in the standard error calculation
  const arma::mat a_matrix = h_matrix * var_ret_inv_cov_ret_fac.t();

  // Center returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Compute intermediate terms for standard error calculation
  const arma::mat fac_cen_gkr_sdf_coeff = factors_centred * gkr_sdf_coefficients;
  const arma::mat ret_cen_var_ret_inv = returns_centred * var_ret_inv;
  const arma::mat ret_cen_fac_cen_gkr_sdf_coeff = returns_centred.each_col() % fac_cen_gkr_sdf_coeff;
  const arma::mat gkr_err = ret_cen_fac_cen_gkr_sdf_coeff.each_row() - mean_returns.t();
  const arma::mat mimicking_err_h =
    (factors_centred - ret_cen_var_ret_inv * covariance_returns_factors) *
    h_matrix;
  const arma::vec u = arma::sum(gkr_err % ret_cen_var_ret_inv, 1);

  // Compute terms used in the standard error calculation
  const arma::mat term1 = gkr_err * a_matrix.t();
  const arma::mat term2 = mimicking_err_h.each_col() % u;

  arma::mat series = term1 + term2;

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(series, hac_prewhite) / std::sqrt(
      static_cast<double>(returns.n_rows)
  );

}
