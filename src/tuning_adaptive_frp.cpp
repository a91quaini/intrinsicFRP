// Author: Alberto Quaini

#include "tuning_adaptive_frp.h"

//////////////////////////////////////////////////////////////
///// MisspecificationIdentificationScoreAdaptiveFRPCpp /////

arma::vec MisspecificationIdentificationScoreAdaptiveFRPCpp(
  const arma::mat& aifrp,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  const double score_no_model = arma::datum::inf;

  // no model will produce a score of infinity
  arma::vec model_score(aifrp.n_cols, arma::fill::value(arma::datum::inf));

  for (unsigned int par = 0; par < aifrp.n_cols; ++par) {

    const arma::uvec idx_selected = arma::find(aifrp.col(par));

    if (! idx_selected.n_elem) continue;

    const arma::mat beta_selected = arma::solve(
      arma::cov(factors),
      covariance_factors_returns.rows(idx_selected),
      arma::solve_opts::likely_sympd
    );

    model_score(par) = HJDistanceCpp(
      beta_selected,
      variance_returns,
      mean_returns
    ) / MinimumSingularValueCpp(beta_selected);

  }

  return model_score;

}
