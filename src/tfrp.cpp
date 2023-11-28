// Author: Alberto Quaini

#include "tfrp.h"
#include "hac_covariance.h"
#include "utils.h"

//////////////////
///// TFRPCpp ////

Rcpp::List TFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  // Check if standard errors should be included in the output
  if (include_standard_errors) {

    // Calculate the covariance between factors and returns
    const arma::mat covariance_factors_returns = arma::cov(factors, returns);
    // Calculate the variance of the returns
    const arma::mat variance_returns = arma::cov(returns);
    // Calculate the mean of the returns
    const arma::vec mean_returns = arma::mean(returns).t();

    // Create a list containing both risk premia and standard errors
    return Rcpp::List::create(
      Rcpp::Named("risk_premia") = TFRPCpp(
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("standard_errors") = StandardErrorsTFRPCpp(
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        hac_prewhite
      )
    );

  } else {

    // If standard errors are not included, only return the risk premia
    return Rcpp::List::create(
      Rcpp::Named("risk_premia") = TFRPCpp(
        arma::cov(factors, returns),
        arma::cov(returns),
        arma::mean(returns).t()
      )
    );
  }

}

arma::vec TFRPCpp(
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  // Return the risk premia computed using the formula Cov(F, R) * Var(R)^(-1) * E(R)
  return covariance_factors_returns * SolveSympd(
    variance_returns,
    mean_returns
  );

}

/////////////////////////////////
///// StandardErrorsTTFRPCpp ////

arma::vec StandardErrorsTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Compute the inverse of variance returns applied to two different matrices
  const arma::mat var_ret_inv = InvSympd(variance_returns);
  const arma::mat var_ret_inv_cov_ret_fac = var_ret_inv *
    covariance_factors_returns.t();

  // Center the returns and factors by subtracting their means
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Compute intermediate terms used in the standard error computation
  const arma::vec ret_cen_var_ret_inv_mean_ret = returns_centred *
    var_ret_inv * mean_returns;
  const arma::mat ret_cen_var_ret_inv_cov_ret_fac = returns_centred *
    var_ret_inv_cov_ret_fac;

  arma::mat series =
    // covariance term
    factors_centred.each_col() % ret_cen_var_ret_inv_mean_ret -
    // variance term
    ret_cen_var_ret_inv_cov_ret_fac.each_col() %
    ret_cen_var_ret_inv_mean_ret +
    // mean term
    returns_centred * var_ret_inv_cov_ret_fac;

  // Calculate the standard errors using HACStandardErrorsCpp function
  return HACStandardErrorsCpp(series, hac_prewhite) /
    std::sqrt(static_cast<double>(returns.n_rows));

}
