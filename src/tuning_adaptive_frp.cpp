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

  // no model will produce a score of infinity
  arma::vec model_score(aifrp.n_cols, arma::fill::value(arma::datum::inf));
  arma::uvec idx_selected_lag;

  for (unsigned int par = 0; par < aifrp.n_cols; ++par) {

    const arma::uvec idx_selected = arma::find(aifrp.col(par));

    // if there are no selected factor, move to the next iteration
    if (! idx_selected.n_elem) continue;

    // if current selected factors are the same as the ones of the
    // previous iteration, set the current model score equal to the previous
    // one and move to the next iteration
    if (idx_selected.n_elem == idx_selected_lag.n_elem) {
      if (arma::all(idx_selected == idx_selected_lag)) {

        model_score(par) = model_score(par - 1);
        continue;

      }
    }

    const arma::mat beta_selected = arma::solve(
      arma::cov(factors.cols(idx_selected)),
      covariance_factors_returns.rows(idx_selected),
      arma::solve_opts::likely_sympd
    ).t();

    model_score(par) = HJDistanceCpp(
      beta_selected,
      variance_returns,
      mean_returns
    ) / MinimumSingularValueCpp(beta_selected);

    idx_selected_lag = idx_selected;

  }

  return model_score;

}


