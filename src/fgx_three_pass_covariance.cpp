// Author: Alberto Quaini

#include "fgx_three_pass_covariance.h"
#include "utils.h"

///////////////////////////////////////
//// CovarianceThreePass //////////////

arma::mat FGXThreePassCovarianceCpp(
  const arma::mat& returns,
  const arma::mat& control_factors,
  const arma::mat& new_factors,
  const arma::vec& sdf_coefficients,
  const arma::uvec& idx_selected
){

  // set number of observations, and lags that are relevant for the Newey-West type
  // covariance estimator
  const unsigned int n_observations = returns.n_rows;
  const unsigned int n_lags = n_observations > 5 ?
    std::floor(
      4 * std::pow(static_cast<double>(0.01 * n_observations), 2.0 / 9.0)
    ) : 0;

  // compute the error matrix in the OLS regressions of new factors on the selected
  // control factors
  const arma::mat selected_controls = control_factors.cols(idx_selected);
  const arma::mat errors = new_factors - selected_controls * SolveSympd(
    selected_controls.t() * selected_controls,
    selected_controls.t() * new_factors
  );

  // compute the errors' precision matrix
  const arma::mat cov_err_inv = InvSympd(errors.t() * errors / n_observations);

  // compute the scaled errors
  const arma::vec scaling = arma::ones(n_observations, 1) -
    arma::join_horiz(new_factors, control_factors) * sdf_coefficients;
  const arma::mat scaled_errors = errors.each_col() % scaling;

  // compute the term of the covariance due to same lag covariances
  const arma::mat first_term = cov_err_inv * scaled_errors.t() * scaled_errors *
    cov_err_inv / static_cast<double>(n_observations);

  // compute the term of the covariance due to lagged covariances
  arma::mat second_term = arma::zeros(new_factors.n_cols, new_factors.n_cols);

  for (int lag = 1; lag <= n_observations; ++lag) {

    const double lag_weight = 1. - static_cast<double>(lag) / (n_lags + 1.0);

    for (int tt = lag; tt < n_observations; ++tt) {

      second_term += lag_weight * (
        scaled_errors.row(tt).t() * scaled_errors.row(tt - lag) +
        scaled_errors.row(tt - lag).t() * scaled_errors.row(tt)
      ) / static_cast<double>(n_observations);

    }

  }

  second_term = cov_err_inv * second_term.t() * second_term * cov_err_inv;

  // return the covariance estimate
  return first_term + second_term;

}
