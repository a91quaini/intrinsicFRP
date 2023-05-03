// Author: Alberto Quaini

#include "frp.h"

/////////////////
///// FRPCpp ////

arma::vec FRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns
) {

  return arma::solve(
    beta.t() * beta, beta.t(), arma::solve_opts::likely_sympd
  ) * mean_returns;

}

////////////////////
///// KRSFRPCpp ////

arma::vec KRSFRPCpp(
  const arma::mat& beta,
  const arma::vec& mean_returns,
  const arma::mat& weighting_matrix
) {

  const arma::mat beta_t_wei_mat_inv = arma::solve(
    weighting_matrix, beta, arma::solve_opts::likely_sympd
  ).t();

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
  const arma::mat& variance_factors,
  const arma::vec& mean_returns,
  const arma::vec& mean_factors
) {

  const arma::mat h_matrix  = arma::inv_sympd(beta.t() * beta);

  const arma::mat a_matrix = h_matrix * beta.t();

  const arma::mat returns_centred = returns.each_row() - mean_returns.t();

  const arma::mat factors_centred = factors.each_row() - mean_factors.t();

  const arma::mat gamma = returns_centred * a_matrix.t();

  const arma::vec gamma_true = a_matrix * mean_returns;

  const arma::mat phi = gamma - factors;

  const arma::mat phi_centred = phi.each_row() - (gamma_true - mean_factors).t();

  const arma::mat fac_centred_var_fac_inv = arma::solve(
    variance_factors,
    factors_centred.t()
  ).t();

  const arma::mat z = fac_centred_var_fac_inv.each_col() % (
    returns_centred * (mean_returns - beta * frp)
  );

  return HACStandardErrorsCpp(
    // mean term
    (gamma.each_row() - (a_matrix * mean_returns).t()) -
    // beta term
    phi_centred.each_col() % (factors_centred * arma::solve(
        variance_factors, gamma_true, arma::solve_opts::likely_sympd
    )) +
    // error term
    (fac_centred_var_fac_inv.each_col() % (
        returns_centred * (mean_returns - beta * frp)
    )) * h_matrix
  );

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
  const arma::mat& variance_factors,
  const arma::vec& mean_returns,
  const arma::vec& mean_factors
) {

  const arma::mat var_ret_inv_beta = arma::solve(
    variance_returns,
    beta,
    arma::solve_opts::likely_sympd
  );

  const arma::mat a_matrix = arma::solve(
    beta.t() * var_ret_inv_beta,
    var_ret_inv_beta.t(),
    arma::solve_opts::likely_sympd
  );

  const arma::mat returns_centred = returns.each_row() - mean_returns.t();

  const arma::mat factors_centred = factors.each_row() - mean_factors.t();

  const arma::mat term1 = returns_centred * a_matrix.t();

  const arma::vec var_ret_inv_mean_ret = arma::solve(
    variance_returns,
    mean_returns,
    arma::solve_opts::likely_sympd
  );

  const arma::mat var_fac_inv = arma::inv_sympd(variance_factors);

  const arma::mat hkrs_var_fac_inv = arma::solve(
    beta.t() * var_ret_inv_beta,
    var_fac_inv,
    arma::solve_opts::likely_sympd
  );

  const arma::vec var_ret_inv_err_krs =
    var_ret_inv_mean_ret - var_ret_inv_beta * krs_frp;

  const arma::mat fac_cen_hkrs_var_fac_inv =
    factors_centred * hkrs_var_fac_inv.t();

  const arma::mat term2 = fac_cen_hkrs_var_fac_inv.each_col() %
    (returns_centred * var_ret_inv_err_krs);

  const arma::mat term3 = term1.each_col() % (returns_centred * var_ret_inv_err_krs);

  const arma::mat akrs_ret_cen_minus_fac_cen = term1 - factors_centred;

  const arma::mat term4 = akrs_ret_cen_minus_fac_cen.each_col() %
    (factors_centred * var_fac_inv * a_matrix * mean_returns);

  return HACStandardErrorsCpp(term1 + term2 - term3 - term4);

}
