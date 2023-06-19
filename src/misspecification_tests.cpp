// Author Alberto Quaini

#include "misspecification_tests.h"

arma::vec2 HJMisspecificationStatisticAndPvalue(
    const arma::mat& returns,
    const arma::mat& factors,
    const arma::mat& beta,
    const arma::mat& variance_returns,
    const arma::vec& mean_returns
) {

  const arma::mat variance_returns_inverse_beta = arma::solve(
    variance_returns, beta, arma::solve_opts::likely_sympd
  );

  const arma::vec krs_rp = arma::solve(
    beta.t() * variance_returns_inverse_beta,
    variance_returns_inverse_beta.t(),
    arma::solve_opts::likely_sympd
  ) * mean_returns;

  const arma::vec pricing_error = mean_returns - beta * krs_rp;

  const arma::vec variance_returns_inverse_pricing_error = arma::solve(
    variance_returns, pricing_error, arma::solve_opts::likely_sympd
  );

  const double hj_distance = arma::dot(
    pricing_error,
    variance_returns_inverse_pricing_error
  );

  const arma::vec u = returns * variance_returns_inverse_pricing_error;

  const arma::vec y =  1. - factors * arma::solve(
    arma::cov(factors, 1), krs_rp, arma::solve_opts::likely_sympd
  );

  const arma::vec q = 2. * u % y - arma::square(u) + hj_distance;

  const double variance_q = ComputeScalarVarianceNeweyWest(q);

  arma::vec2 results;

  results(0) = returns.n_rows * hj_distance * hj_distance / variance_q;

  results(1) = R::pchisq(results(0), 1, false, false);

  return results;

}


