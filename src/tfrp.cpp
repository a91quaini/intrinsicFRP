// Author: Alberto Quaini

#include "tfrp.h"
#include "hac_standard_errors.h"

//////////////////
///// TFRPCpp ////

Rcpp::List TFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool include_standard_errors
) {

  if (include_standard_errors) {

    const arma::mat covariance_factors_returns = arma::cov(factors, returns);
    const arma::mat variance_returns = arma::cov(returns);
    const arma::vec mean_returns = arma::mean(returns).t();

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
        mean_returns
      )
    );

  } else {

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

  return covariance_factors_returns * arma::solve(
    variance_returns,
    mean_returns,
    arma::solve_opts::likely_sympd
  );

}

/////////////////////////////////
///// StandardErrorsTTFRPCpp ////

arma::vec StandardErrorsTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  const arma::mat var_ret_inv_cov_ret_fac = arma::solve(
    variance_returns,
    covariance_factors_returns.t(),
    arma::solve_opts::likely_sympd
  );

  const arma::mat returns_centred = returns.each_row() - mean_returns.t();

  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  const arma::vec ret_cen_var_ret_inv_mean_ret = returns_centred * arma::solve(
    variance_returns,
    mean_returns,
    arma::solve_opts::likely_sympd
  );

  const arma::mat ret_cen_var_ret_inv_cov_ret_fac = returns_centred *
    var_ret_inv_cov_ret_fac;

  return HACStandardErrorsCpp(
    // covariance term
    factors_centred.each_col() % ret_cen_var_ret_inv_mean_ret -
    // variance term
    ret_cen_var_ret_inv_cov_ret_fac.each_col() %
    ret_cen_var_ret_inv_mean_ret +
    // mean term
    returns_centred * var_ret_inv_cov_ret_fac
  );

}
