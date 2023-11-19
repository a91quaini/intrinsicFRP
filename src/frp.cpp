// Author: Alberto Quaini

#include "frp.h"
#include "hac_covariance.h"

/////////////////
///// FRPCpp ////

Rcpp::List FRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust,
  const bool include_standard_errors
) {

  // Check if standard errors are to be included in the calculation
  if (include_standard_errors) {

    // Compute covariance matrices and mean returns
    const arma::mat covariance_factors_returns = arma::cov(factors, returns);
    const arma::mat variance_returns = arma::cov(returns);
    const arma::vec mean_returns = arma::mean(returns).t();

    // Solve for beta coefficients
    const arma::mat beta = arma::solve(
      arma::cov(factors),
      covariance_factors_returns,
      arma::solve_opts::likely_sympd
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
      covariance_factors_returns,
      variance_returns,
      mean_returns
    ) :
      StandardErrorsFRPCpp(
        frp,
        returns,
        factors,
        beta,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      )
    );

  } else {

    // Solve for beta coefficients when standard errors are not included
    const arma::mat beta = arma::solve(
      arma::cov(factors),
      arma::cov(factors, returns),
      arma::solve_opts::likely_sympd
    ).t();

    // Return risk premia only in a list
    return Rcpp::List::create(
      Rcpp::Named("risk_premia") = misspecification_robust ?
    KRSFRPCpp(beta, arma::mean(returns).t(), arma::cov(returns)) :
      FMFRPCpp(beta, arma::mean(returns).t())
    );

  }

}

arma::vec FMFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns
) {

  // Solve for risk premia using Fama-MacBeth method
  return arma::solve(
    beta.t() * beta,
    beta.t(),
    arma::solve_opts::likely_sympd) * mean_returns;

}

////////////////////
///// KRSFRPCpp ////

arma::vec KRSFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns,
  const arma::mat& weighting_matrix
) {

  // Solve for beta coefficients with the weighting matrix
  const arma::mat beta_t_wei_mat_inv = arma::solve(
    weighting_matrix,
    beta,
    arma::solve_opts::likely_sympd
  ).t();

  // Compute risk premia using Kan-Robotti-Shanken method
  return arma::solve(
    beta_t_wei_mat_inv * beta,
    beta_t_wei_mat_inv,
    arma::solve_opts::likely_sympd
  ) * mean_returns;

}

///////////////////////////////
///// StandardErrorsFRPCpp ////

arma::vec StandardErrorsFRPCpp(
  const arma::vec& frp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  // Compute the matrix H as the inverse of beta' * beta
  const arma::mat h_matrix  = arma::inv_sympd(beta.t() * beta);

  // Compute matrix A used in the standard error calculation
  const arma::mat a_matrix = h_matrix * beta.t();

  // Center the returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::vec mean_factors = arma::mean(factors).t();
  const arma::mat factors_centred = factors.each_row() - mean_factors.t();

  // Intermediate matrices for standard error calculation
  const arma::mat gamma = returns_centred * a_matrix.t();
  const arma::vec gamma_true = a_matrix * mean_returns;
  const arma::mat phi = gamma - factors;
  const arma::mat phi_centred = phi.each_row() - (gamma_true - mean_factors).t();

  // Compute variance of factors and solve for fac_centred_var_fac_inv
  const arma::mat variance_factors = arma::cov(factors);
  const arma::mat fac_centred_var_fac_inv = arma::solve(
    variance_factors,
    factors_centred.t(),
    arma::solve_opts::likely_sympd
  ).t();

  // Compute matrix Z used in standard error calculation
  const arma::mat z = fac_centred_var_fac_inv.each_col() % (
    returns_centred * (mean_returns - beta * frp)
  );

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(
    // mean term
    (gamma.each_row() - (a_matrix * mean_returns).t()) -
    // beta term
    phi_centred.each_col() % (factors_centred * arma::solve(
      variance_factors,
      gamma_true,
      arma::solve_opts::likely_sympd
    )) +
    // error term
    (fac_centred_var_fac_inv.each_col() % (
      returns_centred * (mean_returns - beta * frp)
    )) * h_matrix
  ) / std::sqrt(static_cast<double>(returns.n_rows));

}

//////////////////////////////////
///// StandardErrorsKRSFRPCpp ////

arma::vec StandardErrorsKRSFRPCpp(
  const arma::vec& krs_frp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  // Compute the inverse of variance_returns and useful matrices that use it
  const arma::mat var_ret_inv = arma::inv_sympd(variance_returns);
  const arma::mat var_ret_inv_beta = var_ret_inv * beta;
  const arma::vec var_ret_inv_mean_ret = var_ret_inv * mean_returns;

  // Compute matrix A used in the standard error calculation
  const arma::mat a_matrix = arma::solve(
    beta.t() * var_ret_inv_beta,
    var_ret_inv_beta.t(),
    arma::solve_opts::likely_sympd
  );

  // Center returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Compute intermediate terms for standard error calculation
  const arma::vec var_ret_inv_err_krs = var_ret_inv_mean_ret -
    var_ret_inv_beta * krs_frp;
  const arma::mat var_fac_inv = arma::inv_sympd(arma::cov(factors));
  const arma::mat hkrs_var_fac_inv = arma::solve(
    beta.t() * var_ret_inv_beta,
    var_fac_inv,
    arma::solve_opts::likely_sympd
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

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(term1 + term2 - term3 - term4) /
    std::sqrt(static_cast<double>(returns.n_rows));

}
