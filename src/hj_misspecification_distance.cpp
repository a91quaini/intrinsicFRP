// Author Alberto Quaini

#include "hj_misspecification_distance.h"
#include "hac_covariance.h"
#include "utils.h"

///////////////////////////////////////
//// HJMisspecificationDistanceCpp ////

Rcpp::List HJMisspecificationDistanceCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double ci_coverage,
  const bool hac_prewhite
) {

  // Calculate the HJ misspecification statistic and p-value.
  return HJMisspecificationDistanceCpp(
    returns,
    factors,
    arma::cov(returns),
    arma::mean(returns).t(),
    ci_coverage,
    hac_prewhite
  );

}

///////////////////////////////////////
//// HJMisspecificationDistanceCpp ////

Rcpp::List HJMisspecificationDistanceCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double ci_coverage,
  const bool hac_prewhite
) {

  // Compute the inverse variance of returns.
  const arma::mat inv_var_ret = InvSympd(variance_returns);

  // Calculate required vectors and matrices to compute the KRS SDF coefficients.
  const arma::vec var_ret_inv_mean_ret = inv_var_ret * mean_returns;
  const arma::mat covariance_factors_returns = arma::cov(factors, returns);
  const arma::mat var_ret_inv_cov_ret_fac = inv_var_ret *
    covariance_factors_returns.t();

  // Compute the KRS SDF coefficients.
  const arma::vec krs_sdf_coefficients = SolveSympd(
    covariance_factors_returns * var_ret_inv_cov_ret_fac,
    covariance_factors_returns
  ) * var_ret_inv_mean_ret;

  // Calculate the squared Hansen-Jagannathan (HJ) distance.
  const double squared_distance =
    arma::dot(mean_returns, var_ret_inv_mean_ret) - arma::dot(
      mean_returns.t() * var_ret_inv_cov_ret_fac,
      krs_sdf_coefficients
  );

  // Center the factors and calculate the residuals.
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Follow computation in eq. (61) Kan-Robotti (2008) <10.1016/j.jempfin.2008.03.003>.
  const arma::vec u = returns_centred * (
    var_ret_inv_mean_ret - var_ret_inv_cov_ret_fac * krs_sdf_coefficients
  );
  const arma::vec y =  1. - factors_centred * krs_sdf_coefficients;
  arma::vec q = 2. * u % y - arma::square(u) + squared_distance;

  // Compute the HAC standard errors of q.
  const double standard_errors_q = std::sqrt(HACVarianceCpp(q, hac_prewhite));

  // Set the right quantile of the normal based on the ci_coverage.
  const double right_quantile =
    R::qnorm5((1. - ci_coverage) / 2., 0., 1., false, false);

  // Compute the shift around squared_distance determining the confidence interval.
  const double shift = right_quantile * standard_errors_q /
    std::sqrt(static_cast<double>(returns.n_rows));

  // Return the test statistic and the corresponding p-value.
  return Rcpp::List::create(
    Rcpp::Named("squared_distance") = squared_distance,
    Rcpp::Named("lower_bound") = squared_distance - shift,
    Rcpp::Named("upper_bound") = squared_distance + shift
  );

}
