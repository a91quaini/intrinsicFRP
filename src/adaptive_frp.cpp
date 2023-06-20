// Author: Alberto Quaini

#include "adaptive_frp.h"
#include "frp.h"

/////////////////////////////////
///// OptimalAdaptiveFRPCpp /////

Rcpp::List OptimalAdaptiveFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type
) {

  arma::mat aifrp = AdaptiveIFRPCpp(
    IFRPCpp(
      covariance_factors_returns,
      variance_returns,
      mean_returns
    ),
    AdaptiveWeightsCpp(
      returns,
      factors,
      weighting_type
    ),
    penalty_parameters
  );

  const arma::vec model_score = MisspecificationIdentificationScoreAdaptiveFRPCpp(
    aifrp,
    factors,
    covariance_factors_returns,
    variance_returns,
    mean_returns
  );

  unsigned int idx_optimal_parameter = model_score.index_min();

  const arma::uvec idx_selected = arma::find(aifrp.col(idx_optimal_parameter));

  const arma::mat beta_selected = arma::solve(
    arma::cov(factors.cols(idx_selected)),
    covariance_factors_returns.rows(idx_selected),
    arma::solve_opts::likely_sympd
  ).t();

  arma::vec krs_frp(factors.n_cols);
  krs_frp(idx_selected) = KRSFRPCpp(
    beta_selected,
    mean_returns,
    variance_returns
  );

  return Rcpp::List::create(
    Rcpp::Named("risk_premia") = krs_frp,
    Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
    Rcpp::Named("model_score") = model_score
  );

}

/////////////////////////////////////////
///// StandardErrorsAdaptiveFRPCpp /////

arma::vec StandardErrorsAdaptiveFRPCpp(
  const arma::vec& afrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  const arma::uvec idx_selected = arma::find(afrp);

  const arma::mat factors_selected = factors.cols(idx_selected);

  const arma::mat covariance_factors_selected_returns =
    arma::cov(factors_selected, returns);

  arma::vec standard_errors(factors.n_cols);

  standard_errors(idx_selected) = StandardErrorsKRSFRPCpp(
    afrp(idx_selected),
    returns,
    factors_selected,
    arma::solve(
      arma::cov(factors_selected),
      covariance_factors_selected_returns,
      arma::solve_opts::likely_sympd
    ).t(),
    covariance_factors_selected_returns,
    variance_returns,
    arma::cov(factors_selected),
    mean_returns,
    arma::mean(factors_selected).t()
  );

  return standard_errors;

}
