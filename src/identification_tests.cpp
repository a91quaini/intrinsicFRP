// Author: Alberto Quaini

#include "identification_tests.h"
#include "utils.h"

namespace {

arma::mat SolveSquareSystem(
  const arma::mat& A,
  const arma::mat& b
) {

  try {
    return arma::solve(A, b, arma::solve_opts::no_approx);
  } catch (const std::runtime_error&) {
    return arma::solve(A, b, arma::solve_opts::force_approx);
  }

}

} // namespace

///////////////////////////////////////
///// ChenFang2019BetaRankTestCpp /////

Rcpp::List ChenFang2019BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const unsigned int n_bootstrap,
  const double target_level_kp2006_rank_test
) {

  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;

  if (n_factors >= n_returns) {
    Rcpp::stop("Chen Fang 2019 test: n_factors must be < n_returns");
  }

  const unsigned int n_observations = returns.n_rows;
  const arma::mat beta = ScaledFactorLoadingsCpp(returns, factors).t();

  arma::mat U(n_returns, n_returns);
  arma::mat V(n_factors, n_factors);
  arma::vec sv(n_factors);
  arma::svd(U, sv, V, beta);

  const double statistic = static_cast<double>(n_observations) *
    std::pow(sv(n_factors - 1), 2.0);

  const unsigned int rank_estimate = target_level_kp2006_rank_test > 0. ?
    Rcpp::as<unsigned int>(
      IterativeKleibergenPaap2006BetaRankTestCpp(
        returns,
        factors,
        target_level_kp2006_rank_test
      )["rank"]
    ) :
    arma::sum(sv >= std::pow(static_cast<double>(n_observations), -1.0 / 4.0));

  if (rank_estimate == n_factors) {
    return Rcpp::List::create(
      Rcpp::Named("statistic") = statistic,
      Rcpp::Named("p-value") = 0.
    );
  }

  const arma::mat U22 = U.tail_cols(n_returns - rank_estimate);
  const arma::mat V22 = V.tail_cols(n_factors - rank_estimate);
  arma::vec min_sv_boot(n_bootstrap);

  const double sqrt_n_obs = std::sqrt(static_cast<double>(n_observations));

  for (unsigned int boot = 0; boot < n_bootstrap; ++boot) {

    const arma::uvec boot_indices = arma::randi<arma::uvec>(
      n_observations,
      arma::distr_param(0, n_observations - 1)
    );

    const arma::mat beta_boot = ScaledFactorLoadingsCpp(
      returns.rows(boot_indices),
      factors.rows(boot_indices)
    ).t();

    min_sv_boot(boot) = arma::min(arma::svd(
      sqrt_n_obs * U22.t() * (beta_boot - beta) * V22
    ));

  }

  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("p-value") = static_cast<double>(arma::sum(
      min_sv_boot % min_sv_boot >= statistic
    )) / static_cast<double>(n_bootstrap)
  );

}


///////////////////////////////////////////////////////////////
///// KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp /////

arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::vec& theta_vectorised,
  const arma::mat& U,
  const arma::mat& V,
  const arma::mat& W,
  const unsigned int q
) {

  const unsigned int n_factors = U.n_rows;
  const unsigned int n_returns = V.n_cols;

  const arma::mat U22 = U.submat(q, q, n_factors - 1, n_factors - 1);
  const arma::mat V22 = V.submat(q, q, n_returns - 1, n_returns - 1);

  arma::vec U22_eval;
  arma::vec V22_eval;
  arma::mat U22_evec;
  arma::mat V22_evec;

  arma::eig_sym(U22_eval, U22_evec, U22 * U22.t());
  arma::eig_sym(V22_eval, V22_evec, V22 * V22.t());

  const arma::vec sqrt_U22_eval = arma::sqrt(
    arma::clamp(U22_eval, 0.0, arma::datum::inf)
  );
  const arma::vec sqrt_V22_eval = arma::sqrt(
    arma::clamp(V22_eval, 0.0, arma::datum::inf)
  );

  const arma::mat sqrt_U22 = U22_evec * arma::diagmat(sqrt_U22_eval) * U22_evec.t();
  const arma::mat sqrt_V22 = V22_evec * arma::diagmat(sqrt_V22_eval) * V22_evec.t();

  const arma::mat A_qperp = U.tail_cols(n_factors - q) *
    SolveSquareSystem(U22, sqrt_U22);
  const arma::mat B_qperp = sqrt_V22 *
    SolveSquareSystem(V22.t(), V.tail_cols(n_returns - q).t());

  const arma::mat kron_BA_qperp = arma::kron(B_qperp, A_qperp.t());
  const arma::vec lambda_q = kron_BA_qperp * theta_vectorised;
  const arma::mat Omega_q = kron_BA_qperp * W * kron_BA_qperp.t();

  const double statistic = arma::dot(lambda_q, SolveSympd(Omega_q, lambda_q));
  const double pvalue = R::pchisq(
    statistic,
    static_cast<double>((n_factors - q) * (n_returns - q)),
    false,
    false
  );

  return arma::vec2{statistic, pvalue};

}

//////////////////////////////////////////////////////
///// IterativeKleibergenPaap2006BetaRankTestCpp /////

Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level
) {

  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;
  const unsigned int n_observations = returns.n_rows;

  if (n_factors >= n_returns) {
    Rcpp::stop("Kleibergen Paap test: n_factors must be < n_returns");
  }

  const arma::mat fac_t_fac = factors.t() * factors;
  const arma::mat pi_ols = SolveSympd(fac_t_fac, factors.t() * returns);
  const arma::mat theta = ScaledFactorLoadingsCpp(returns, factors);

  const arma::mat U_ret = arma::chol(returns.t() * returns);
  const arma::mat G_fac = arma::chol(fac_t_fac);

  const arma::mat residuals = returns - factors * pi_ols - arma::repmat(
    arma::mean(returns),
    n_observations,
    1
  );
  const arma::mat err1 = arma::solve(
    arma::trimatl(U_ret.t()),
    residuals.t()
  ).t();
  const arma::mat err2 = arma::solve(
    arma::trimatl(G_fac.t()),
    factors.t()
  ).t();
  const arma::mat err = arma::repelem(err1, 1, n_factors) %
    arma::repmat(err2, 1, n_returns);
  const arma::mat W = err.t() * err;

  arma::mat U(n_factors, n_factors);
  arma::mat V(n_returns, n_returns);
  arma::vec sv(n_factors);
  arma::svd(U, sv, V, theta);

  arma::mat output(2, n_factors, arma::fill::zeros);

  for (unsigned int q = 0; q < n_factors; ++q) {
    output.col(q) = KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
      arma::vectorise(theta),
      U,
      V,
      W,
      q
    );
  }

  const arma::uvec idx_accept = arma::find(
    output.row(1) > target_level / static_cast<double>(n_factors)
  );
  const unsigned int rank = idx_accept.empty() ? n_factors : arma::min(idx_accept);

  return Rcpp::List::create(
    Rcpp::Named("rank") = rank,
    Rcpp::Named("q") = arma::regspace<arma::rowvec>(0, n_factors - 1),
    Rcpp::Named("statistics") = output.row(0),
    Rcpp::Named("pvalues") = output.row(1)
  );

}

///////////////////////////////////
///// ScaledFactorLoadingsCpp /////

arma::mat ScaledFactorLoadingsCpp(
  const arma::mat& returns,
  const arma::mat& factors
) {

  const arma::mat L_fac = arma::chol(factors.t() * factors, "lower");
  const arma::mat U_ret = arma::chol(returns.t() * returns);
  const arma::mat pi = arma::solve(
    arma::trimatl(L_fac),
    factors.t() * returns
  );

  return pi * arma::inv(arma::trimatu(U_ret));

}
