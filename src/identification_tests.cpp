// Author: Alberto Quaini

#include "identification_tests.h"

double MinSingularValue(const arma::mat matrix){

  arma::vec sv = arma::svd(matrix);

  return sv(sv.n_elem - 1);

}

arma::vec2 BetaRankChenFang2019StatisticAndPvalue(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& beta,
  const int sv_threshold_type,
  const unsigned int n_bootstrap
) {

  const unsigned int n_observations = returns.n_rows;
  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;
  const unsigned int n_sv = std::min(n_returns, n_factors);

  // SVD of beta
  arma::mat left(n_returns, n_returns);
  arma::mat right(n_factors, n_factors);
  arma::vec sv(n_sv);

  arma::svd(left, sv, right, beta);

  arma::vec2 output;
  output(0) = n_observations * sv(n_sv - 1) * sv(n_sv - 1);

  const double sv_threshold = sv_threshold_type == 0 ?
    std::pow(n_observations, -1./4.) :
    std::pow(n_observations, -1./3.);

  const unsigned int rank_estimate = arma::sum(sv >= sv_threshold);

  if (rank_estimate == n_sv) {
    output(1) = 0.;
    return output;
  }

  const arma::mat left_zero = left.tail_cols(n_returns - rank_estimate);
  const arma::mat right_zero = right.tail_cols(n_factors - rank_estimate);
  arma::vec sv_boot(n_bootstrap);

  for (unsigned int boot = 0; boot < n_bootstrap; ++boot) {

    const arma::uvec boot_indexes = arma::randi<arma::uvec>(
      n_observations,
      arma::distr_param(0, n_observations - 1)
    );

    const arma::mat factors_boot = factors.rows(boot_indexes);

    const arma::mat beta_boot = arma::solve(
      arma::cov(factors_boot, 1),
      arma::cov(factors_boot, returns.rows(boot_indexes), 1),
      arma::solve_opts::likely_sympd
    ).t();

    const arma::vec sv_statistic_boot = arma::svd(
      std::sqrt(n_observations) * left_zero.t() * (beta_boot - beta) * right_zero
    );

    sv_boot(boot) = sv_statistic_boot(sv_statistic_boot.n_elem - 1);

  }

  output(1) = (double)arma::sum(sv_boot % sv_boot >= output(0)) / n_bootstrap;

  return output;

}
