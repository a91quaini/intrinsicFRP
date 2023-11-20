// Author Alberto Quaini

#include "misspecification_tests.h"
#include "hac_covariance.h"

////////////////////////////////////
//// HJMisspecificationTestCpp ////

Rcpp::List HJMisspecificationTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool hac_prewhite
) {

  // Calculate the HJ misspecification statistic and p-value.
  return HJMisspecificationStatisticAndPvalueCpp(
    returns,
    factors,
    arma::cov(returns),
    arma::mean(returns).t(),
    hac_prewhite
  );

}

/////////////////////////////////////////////////
//// HJMisspecificationStatisticAndPvalueCpp ////

Rcpp::List HJMisspecificationStatisticAndPvalueCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Perform singular value decomposition on the variance of returns.
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, variance_returns);

  // Compute different powers of variance_returns.
  const arma::mat inv_var_ret = eigvec * arma::diagmat(1. / eigval) * eigvec.t();

  // Calculate required vectors and matrices to compure the KRS SDF coefficients.
  const arma::vec var_ret_inv_mean_ret = arma::solve(
    variance_returns,
    mean_returns,
    arma::solve_opts::likely_sympd
  );
  const arma::mat covariance_factors_returns = arma::cov(factors, returns);
  const arma::mat var_ret_inv_cov_ret_fac = arma::solve(
    variance_returns,
    covariance_factors_returns.t(),
    arma::solve_opts::likely_sympd
  );

  // Compute the KRS SDF coefficients.
  const arma::vec krs_sdf_coefficients = arma::solve(
    covariance_factors_returns * var_ret_inv_cov_ret_fac,
    covariance_factors_returns,
    arma::solve_opts::likely_sympd
  ) * var_ret_inv_mean_ret;

  // Calculate the Hansen-Jagannathan (HJ) statistic.
  const double hj_distance = arma::dot(mean_returns, var_ret_inv_mean_ret) -
    arma::dot(
      mean_returns.t() * var_ret_inv_cov_ret_fac,
      krs_sdf_coefficients
  );

  // Center the factors and calculate the residuals.
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Follow computation in eq. (61) Kan-Robotti (2008)
  // <10.1016/j.jempfin.2008.03.003>.
  const arma::vec u = returns_centred * (
    var_ret_inv_mean_ret - var_ret_inv_cov_ret_fac * krs_sdf_coefficients
  );
  const arma::vec y =  1. - factors_centred * krs_sdf_coefficients;
  arma::vec q = 2. * u % y - arma::square(u) + hj_distance;

  const double variance_q = HACVarianceCpp(q, hac_prewhite);

  // Compute the squared standardized HJ test statistic.
  const double statistic = returns.n_rows * hj_distance * hj_distance / variance_q;

  // Return the test statistic and the corresponding p-value.
  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("p-value") = R::pchisq(statistic, 1, false, false)
  );

}
