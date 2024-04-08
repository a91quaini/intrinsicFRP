// Author: Alberto Quaini

#include "frp.h"
#include "hac_covariance.h"
#include "gkr_factor_screening.h"
#include "utils.h"

/////////////////
///// FRPCpp ////

Rcpp::List FRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust,
  const bool include_standard_errors,
  const bool hac_prewhite,
  const double target_level_gkr2014_screening
) {

  // If target_level_gkr2014_screening > 0, perform the GKR screening procedure
  if (target_level_gkr2014_screening > 0.) {

    // Store the indices of the selected factors
    const arma::uvec selected_factor_indices = static_cast<arma::uvec>(
      GKRFactorScreeningCpp(
        returns,
        factors,
        target_level_gkr2014_screening,
        hac_prewhite
      )["selected_factor_indices"]
    );

    // If the GKR procedure removed all the factors, return a list containing empty vectors.
    if (selected_factor_indices.empty()) return include_standard_errors ?
      Rcpp::List::create(
        Rcpp::Named("risk_premia") = arma::vec(),
        Rcpp::Named("standard_errors") = arma::vec(),
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      ) :
      Rcpp::List::create(
        Rcpp::Named("risk_premia") = arma::vec(),
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      );

    // Otherwise, compute the factor risk premia and eventual standard errors
    const Rcpp::List output = ReturnFRPCpp(
      returns,
      factors.cols(selected_factor_indices),
      misspecification_robust,
      include_standard_errors,
      hac_prewhite
    );

    // Return the output list with additionally the selected factor indices
    return include_standard_errors ?
      Rcpp::List::create(
        Rcpp::Named("risk_premia") = output["risk_premia"],
        Rcpp::Named("standard_errors") = output["standard_errors"],
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      ) :
      Rcpp::List::create(
        Rcpp::Named("risk_premia") = output["risk_premia"],
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      );

  }

  // Otherwise, return list of results without selected factor indices
  return ReturnFRPCpp(
    returns,
    factors,
    misspecification_robust,
    include_standard_errors,
    hac_prewhite
  );

}

////////////////////////
////// ReturnFRPCpp ////

Rcpp::List ReturnFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust,
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  // Check if standard errors are to be included in the calculation
  if (include_standard_errors) {

    // Compute covariance matrices and mean returns
    const arma::mat covariance_factors_returns = arma::cov(factors, returns);
    const arma::mat variance_returns = arma::cov(returns);
    const arma::vec mean_returns = arma::mean(returns).t();

    // Solve for beta coefficients
    const arma::mat beta = SolveSympd(
      arma::cov(factors),
      covariance_factors_returns
    ).t();

    // Compute risk premia based on whether misspecification is robust or not
    const arma::vec frp = misspecification_robust ?
    KRSFRPCpp(beta, mean_returns, variance_returns) :
      FMFRPCpp(beta, mean_returns);

    // Return risk premia and standard errors in a list
    return Rcpp::List::create(
        Rcpp::Named("risk_premia") =  frp,
        Rcpp::Named("standard_errors") = misspecification_robust ?
          StandardErrorsKRSFRPCpp(
            frp,
            returns,
            factors,
            beta,
            variance_returns,
            mean_returns,
            hac_prewhite
          ) :
          StandardErrorsFRPCpp(
            frp,
            returns,
            factors,
            beta,
            variance_returns,
            mean_returns,
            hac_prewhite
          )
    );

  } else {

    // Solve for beta coefficients when standard errors are not included
    const arma::mat beta = SolveSympd(
      arma::cov(factors),
      arma::cov(factors, returns)
    ).t();

    // Return risk premia only in a list
    return Rcpp::List::create(
      Rcpp::Named("risk_premia") = misspecification_robust ?
        KRSFRPCpp(beta, arma::mean(returns).t(), arma::cov(returns)) :
        FMFRPCpp(beta, arma::mean(returns).t())
    );

  }

}

////////////////////
////// FMFRPCpp ////

arma::vec FMFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns
) {

  // Solve for risk premia using Fama-MacBeth method
  return SolveSympd(beta.t() * beta, beta.t()) * mean_returns;

}

////////////////////
///// KRSFRPCpp ////

arma::vec KRSFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns,
  const arma::mat& weighting_matrix
) {

  // Solve for beta coefficients with the weighting matrix
  const arma::mat beta_t_wei_mat_inv = SolveSympd(
    weighting_matrix,
    beta
  ).t();

  // Compute risk premia using Kan-Robotti-Shanken method
  return SolveSympd(
    beta_t_wei_mat_inv * beta,
    beta_t_wei_mat_inv
  ) * mean_returns;

}

///////////////////////////////
///// StandardErrorsFRPCpp ////

arma::vec StandardErrorsFRPCpp(
  const arma::vec& frp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Compute the matrix H as the inverse of beta' * beta
  const arma::mat h_matrix  = InvSympd(beta.t() * beta);

  // Center the returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // term1: mean returns
  const arma::mat term1 = returns_centred * beta * h_matrix;

  // term2: beta
  const arma::mat phi_centred = term1 - factors_centred;
  const arma::mat fac_cen_var_fac_inv = SolveSympd(
    arma::cov(factors),
    factors_centred.t()
  ).t();

  // term3: misspecification
  const arma::mat z = fac_cen_var_fac_inv.each_col() % (
    returns_centred * (mean_returns - beta * frp)
  );

  arma::mat series =
    // mean term
    term1 -
    // beta term
    phi_centred.each_col() % (fac_cen_var_fac_inv * frp) +
    // misspecification term
    z * h_matrix;

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(series, hac_prewhite) /
    std::sqrt(static_cast<double>(returns.n_rows));

}

//////////////////////////////////
///// StandardErrorsKRSFRPCpp ////

arma::vec StandardErrorsKRSFRPCpp(
  const arma::vec& krs_frp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Compute the inverse of variance_returns and useful matrices that use it
  const arma::mat var_ret_inv = InvSympd(variance_returns);
  const arma::mat var_ret_inv_beta = var_ret_inv * beta;
  const arma::vec var_ret_inv_mean_ret = var_ret_inv * mean_returns;

  // Compute matrix A used in the standard error calculation
  const arma::mat a_matrix = SolveSympd(
    beta.t() * var_ret_inv_beta,
    var_ret_inv_beta.t()
  );

  // Center returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Compute intermediate terms for standard error calculation
  const arma::vec var_ret_inv_err_krs = var_ret_inv_mean_ret -
    var_ret_inv_beta * krs_frp;
  const arma::mat var_fac_inv = InvSympd(arma::cov(factors));
  const arma::mat hkrs_var_fac_inv = SolveSympd(
    beta.t() * var_ret_inv_beta,
    var_fac_inv
  );

  // Compute terms used in the standard error calculation
  const arma::mat term1 = returns_centred * a_matrix.t();
  const arma::mat fac_cen_hkrs_var_fac_inv =
    factors_centred * hkrs_var_fac_inv.t();
  const arma::mat term2 = fac_cen_hkrs_var_fac_inv.each_col() %
    (returns_centred * var_ret_inv_err_krs);
  const arma::mat term3 = term1.each_col() % (returns_centred * var_ret_inv_err_krs);
  const arma::mat akrs_ret_cen_minus_fac_cen = term1 - factors_centred;
  const arma::mat term4 = akrs_ret_cen_minus_fac_cen.each_col() %
    (factors_centred * var_fac_inv * a_matrix * mean_returns);

  arma::mat series = term1 + term2 - term3 - term4;

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(series, hac_prewhite) /
    std::sqrt(static_cast<double>(returns.n_rows));

}
