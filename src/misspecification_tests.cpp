// Author Alberto Quaini

#include "misspecification_tests.h"
#include "hac_standard_errors.h"

////////////////////////////////////
//// HJMisspecificationTestCpp ////

Rcpp::List HJMisspecificationTestCpp(
  const arma::mat& returns,
  const arma::mat& factors
) {

  return HJMisspecificationStatisticAndPvalueCpp(
    returns,
    factors,
    arma::solve(
      arma::cov(factors),
      arma::cov(factors, returns),
      arma::solve_opts::likely_sympd
    ).t(),
    arma::cov(returns),
    arma::mean(returns).t()
  );

}

/////////////////////////////////////////////////
//// HJMisspecificationStatisticAndPvalueCpp ////

Rcpp::List HJMisspecificationStatisticAndPvalueCpp(
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

  const double variance_q = NeweyWestVarianceOfScalarSeriesCpp(q);

  const double statistic = returns.n_rows * hj_distance * hj_distance / variance_q;

  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("p-value") = R::pchisq(statistic, 1, false, false)
  );

}

// double HJDistanceCpp(
//   const arma::mat& beta,
//   const arma::mat& variance_returns,
//   const arma::vec& mean_returns
// ) {
//
//   const arma::mat variance_returns_inverse_beta = arma::solve(
//     variance_returns, beta, arma::solve_opts::likely_sympd
//   );
//
//   const arma::vec krs_rp = arma::solve(
//     beta.t() * variance_returns_inverse_beta,
//     variance_returns_inverse_beta.t(),
//     arma::solve_opts::likely_sympd
//   ) * mean_returns;
//
//   const arma::vec pricing_error = mean_returns - beta * krs_rp;
//
//   const arma::vec variance_returns_inverse_pricing_error = arma::solve(
//     variance_returns, pricing_error, arma::solve_opts::likely_sympd
//   );
//
//   return arma::dot(
//     pricing_error,
//     variance_returns_inverse_pricing_error
//   );
//
// }

