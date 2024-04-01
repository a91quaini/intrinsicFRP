// Author: Alberto Quaini

#include "gkr_factor_screening.h"
#include "sdf_coefficients.h"
#include "hac_covariance.h"
#include "frp.h"
#include "utils.h"

////////////////////////////////
///// GKRFactorScreeningCpp ////

Rcpp::List GKRFactorScreeningCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level,
  const bool hac_prewhite
) {

  return GKRFactorScreeningCpp(
    returns,
    factors,
    arma::cov(returns, factors),
    arma::cov(returns),
    arma::mean(returns).t(),
    target_level,
    hac_prewhite
  );

}

////////////////////////////////
///// GKRFactorScreeningCpp ////

Rcpp::List GKRFactorScreeningCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double target_level,
  const bool hac_prewhite
) {

  // Copy matrices holding data that may change after the removal of a factor.
  arma::mat factors_new = factors;
  arma::mat covariance_returns_factors_new = covariance_returns_factors;

  // Store the indices of the selected factors
  arma::uvec selected_factor_indices = arma::regspace<arma::uvec>(0, factors.n_cols - 1);

  // Loop over potential factor removals.
  for (int idx = 0; idx < factors.n_cols; ++idx) {

    // Compute the critical value of the normal distribution as the quantile
    // corresponding to the CDF evaluated at  1 - (target_level / n_factors).
    // Note that the Bonferroni correction adjusts to the current number of factors.
    const double critical_value = R::qchisq(
      1 - (target_level / factors_new.n_cols),
      1.0,
      true,
      false
    );
    // Rcpp::Rcout << "critical value = " << critical_value << "\n";
    //
    // const arma::mat beta = arma::solve(
    //   arma::cov(factors_new), covariance_returns_factors_new.t(),
    //   arma::solve_opts::likely_sympd
    // ).t();
    // Rcpp::Rcout << "KRS FRP = " << KRSFRPCpp(
    //     beta, mean_returns, variance_returns
    // ) << "\n";
    // Rcpp::Rcout << "KRS FRP standard errors = " << StandardErrorsFRPCpp(
    //     KRSFRPCpp(
    //       beta, mean_returns, variance_returns
    //     ), returns, factors_new, beta, covariance_returns_factors_new.t(), variance_returns, mean_returns
    // ) << "\n";

    // Compute the GKR SDF coefficients and their standard errors.
    const arma::vec gkr_sdf_coefficients = GKRSDFCoefficientsCpp(
      covariance_returns_factors_new,
      variance_returns,
      mean_returns
    );
    // Rcpp::Rcout << "SDF coefficients = " << gkr_sdf_coefficients << "\n";
    const arma::vec standard_errors = StandardErrorsGKRSDFCoefficientsCpp(
      gkr_sdf_coefficients,
      returns,
      factors_new,
      covariance_returns_factors_new,
      variance_returns,
      mean_returns,
      hac_prewhite
    );
    // Rcpp::Rcout << "standard_errors = " << standard_errors << "\n";

    // Compute the t-statistic for each factor risk premia.
    const arma::vec t_stat = gkr_sdf_coefficients / standard_errors;
    const arma::vec squared_t_stat = arma::square(t_stat);
    // Rcpp::Rcout << "squared_t_stat = " << squared_t_stat << "\n";

    // Store the index of the minimum absolute t-statistic
    const unsigned int idx_min = arma::index_min(squared_t_stat);

    // If the minimum t-statistic is above the critical value,
    // return the results.
    if (squared_t_stat(idx_min) > critical_value) {

      // Return a list contaning the GKR SDF coefficients, standard errors and
      // t-statistics. Further return the indices of the selected factors.
      return Rcpp::List::create(
        Rcpp::Named("sdf_coefficients") = gkr_sdf_coefficients,
        Rcpp::Named("standard_errors") = standard_errors,
        Rcpp::Named("t_statistics") = t_stat,
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      );

    }

    // Otherwise, remove the factor associated to the minimum t-statistic
    // from the data.
    factors_new.shed_col(idx_min);

    // if n_factors is zero, exit the loop to avoid meaningless computations.
    if (factors_new.empty()) break;
    covariance_returns_factors_new.shed_col(idx_min);
    selected_factor_indices.shed_row(idx_min);

  }

  // If the GKR procedure removed all the factors, return a list containing empty vectors.
  return Rcpp::List::create(
    Rcpp::Named("sdf_coefficients") = arma::vec(),
    Rcpp::Named("standard_errors") = arma::vec(),
    Rcpp::Named("t_statistics") = arma::vec(),
    Rcpp::Named("selected_factor_indices") = arma::uvec()
  );

}
