// Author: Alberto Quaini

#include "rank_tests.h"

Rcpp::List BetaRankKleibergenPaap2006StatisticsAndPvaluesCpp(
  const arma::mat& returns,
  const arma::mat& factors
) {

  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;

  if (n_factors >= n_returns) Rcpp::stop("KP test: n_factors must be < n_returns");

  const unsigned int n_observations = factors.n_rows;

  const arma::mat fac_t_fac = factors.t() * factors;
  const arma::mat chol_ret_t_ret = arma::chol(returns.t() * returns);
  const arma::mat chol_fac_t_fac = arma::chol(fac_t_fac);

  const arma::mat pi = arma::solve(
    fac_t_fac, factors.t() * returns,
    arma::solve_opts::likely_sympd
  );

  const arma::mat theta = chol_fac_t_fac * arma::solve(
    chol_ret_t_ret, pi, arma::solve_opts::likely_sympd
  ).t();

  arma::mat U(n_factors, n_factors);
  arma::mat V(n_returns, n_returns);
  arma::vec s(n_factors);

  arma::svd(U, s, V, theta);

  // White covariance matrix estimator with scaling matrices incorporated
  arma::mat V_hat(n_returns * n_factors, n_returns * n_factors);

  const arma::mat residuals = returns - factors * pi - arma::repmat(
    arma::mean(returns), n_observations, 1
  );
  const arma::mat err1 = arma::solve(
    chol_ret_t_ret, residuals, arma::solve_opts::likely_sympd
  ).t();
  const arma::mat err2 = arma::solve(
    chol_fac_t_fac, factors, arma::solve_opts::likely_sympd
  ).t();

  for (unsigned int obs = 0; obs < n_observations; ++obs) {

    const arma::vec err = arma::kron(err1, err2.row(obs)).t();
    V_hat += err * err.t();

  }

  const arma::vec vec_theta = arma::vectorise(theta);

  arma::vec statistics(n_factors - 1);
  arma::vec pvalues(n_factors - 1);

  // loop
  for (unsigned int rank = 0; rank < n_factors - 2; rank++) {

    const arma::mat U22 = U.submat(rank + 1, n_factors - 1, rank + 1, n_factors - 1);
    const arma::mat V22 = V.submat(rank + 1, n_returns - 1, rank + 1, n_returns - 1);

    const arma::mat sqrt_U22 = arma::sqrtmat_sympd(U22 * U22.t());
    const arma::mat sqrt_V22 = arma::sqrtmat_sympd(V22 * V22.t());

    // equation (12) in the Kleibergen Paap 2006 paper
    const arma::mat A_perp = U.tail_cols(n_factors - rank) * arma::solve(
      U22, sqrt_U22, arma::solve_opts::likely_sympd
    );
    const arma::mat B_perp = sqrt_V22 * arma::solve(
      V22.t(), V.tail_cols(n_factors - rank).t()
    );

    const arma::mat kron_BA_perp = arma::kron(B_perp, A_perp.t());

    // equation below (21) in the Kleibergen Paap 2006 paper
    const arma::vec lambda = kron_BA_perp * vec_theta;

    // covariance matrix
    const arma::mat V_lambda = kron_BA_perp * V_hat * kron_BA_perp.t();

    // test statistic
    statistics(rank) = arma::dot(lambda, arma::solve(
      V_lambda, lambda, arma::solve_opts::likely_sympd
    ));

    // pvalue
    pvalues(rank) = 1. - R::pchisq(
      statistics(rank),
      (n_factors - rank) * (n_returns - rank),
      true,
      false
    );

  }

  return Rcpp::List::create(
    Rcpp::Named("statistics") = statistics,
    Rcpp::Named("pvalues") = pvalues
  );

}

arma::vec2 BetaRankChenFang2019StatisticAndPvalueCpp(
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

  // set singular value threshold: a number that tends to 0 as n_observations
  // grows, but tends to infinity if multiplied by sqrt(n_observations)
  const double sv_threshold = sv_threshold_type == 0 ?
    std::pow(n_observations, -1./3) / (std::log(n_factors) * std::log(n_returns)) :
    std::pow(n_observations, -1./3);

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
