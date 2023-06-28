// Author: Alberto Quaini

#include "rank_tests.h"

///////////////////////////////////////////////////////////////
///// KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp /////

arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::mat& theta,
  const arma::mat& U,
  const arma::mat& V,
  const arma::mat& W,
  const unsigned int q,
  const unsigned int n_observations
) {

  const unsigned int n_factors = theta.n_rows;
  const unsigned int n_returns = theta.n_cols;

  // see eq. (3) in Kleibergen Paap 2006
  const arma::mat U22 = U.submat(q, q, n_factors - 1, n_factors - 1);
  const arma::mat V22 = V.submat(q, q, n_factors - 1, n_factors - 1);

  Rcpp::Rcout << "V22 = " << V22 << "\n";

  const arma::mat sqrt_U22 = arma::sqrtmat_sympd(U22 * U22.t());
  const arma::mat sqrt_V22 = arma::sqrtmat_sympd(V22 * V22.t());

  Rcpp::Rcout << "sqrt_U22 = " << sqrt_U22 << "\n";

  // see eq. (12) in Kleibergen Paap 2006
  const arma::mat A_qperp = U.tail_cols(n_factors - q) * arma::solve(
    U22, sqrt_U22,
    arma::solve_opts::likely_sympd
  );
  const arma::mat B_qperp = sqrt_V22 * arma::solve(
    V22.t(), V.tail_cols(n_factors - q).t(),
    arma::solve_opts::likely_sympd
  );

  Rcpp::Rcout << "A_qperp = " << A_qperp << "\n";
  Rcpp::Rcout << "B_qperp = " << B_qperp << "\n";

  const arma::mat kron_BA_qperp = arma::kron(B_qperp, A_qperp.t());
  Rcpp::Rcout << "kron_BA_qperp = " << kron_BA_qperp.submat(0,2,0,2) << "\n";


  // see below eq. (21) in Kleibergen Paap 2006
  // use that vec(A'lB') = (B kron A') vec(l)
  const arma::vec lambda_q = kron_BA_qperp * arma::vectorise(theta);
  const arma::mat Omega_q = kron_BA_qperp * W * kron_BA_qperp.t();

  Rcpp::Rcout << "lambda_q = " << lambda_q(0,2) << "\n";
  Rcpp::Rcout << "Omega_q = " << Omega_q.submat(0,2,0,2) << "\n";

  arma::vec2 output;

  // test statistic
  output(0) = n_observations * arma::dot(lambda_q, arma::solve(
      Omega_q, lambda_q,
      arma::solve_opts::likely_sympd
  ));


  // p-value
  output(1) = 1. - R::pchisq(
    output(0),
    (n_factors - q) * (n_returns - q),
    true,
    false
  );
  Rcpp::Rcout << "output(0)  = " << output(0)  << "\n";
  Rcpp::Rcout << "output(1)  = " << output(1)  << "\n";

  return output;

}

//////////////////////////////////////////////////////
///// IterativeKleibergenPaap2006BetaRankTestCpp /////

Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double level,
  const bool scaling
) {

  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;

  if (n_factors >= n_returns) Rcpp::stop("Kleibergen Paap test: n_factors must be < n_returns");

  const unsigned int n_observations = factors.n_rows;

  const arma::mat fac_t_fac = factors.t() * factors;

  // see below eq. (16) in Kleibergen Paap 2006
  const arma::mat pi = arma::solve(
    fac_t_fac, factors.t() * returns,
    arma::solve_opts::likely_sympd
  );
  Rcpp::Rcout << "pi = " << pi << "\n";


  // see eq. (16) in Kleibergen Paap 2006
  const arma::mat F_t_inv = arma::chol(returns.t() * returns);
  Rcpp::Rcout << "F_t_inv = " << F_t_inv << "\n";
  const arma::mat G = arma::chol(fac_t_fac);
  const arma::mat theta = G * arma::solve(
    arma::trimatl(F_t_inv.t()), pi.t()
  ).t();
  Rcpp::Rcout << "theta = " << theta << "\n";

  // White covariance matrix estimator with scaling matrices incorporated
  arma::mat W(n_returns * n_factors, n_returns * n_factors, arma::fill::eye);

  const arma::mat residuals = returns - factors * pi;// - arma::repmat(
  //  arma::mean(returns), n_observations, 1
  //);
  Rcpp::Rcout << "residuals = " << residuals.submat(0,0,2,2) << "\n";
  const arma::mat err1 = arma::solve(
    arma::trimatl(F_t_inv.t()), residuals.t()
  ).t();
  Rcpp::Rcout << "err1 = " << err1.submat(0,0,2,2) << "\n";
  const arma::mat err2 = arma::solve(
    arma::trimatl(G.t()), factors.t()
  ).t();
  Rcpp::Rcout << "err2 = " << err2.submat(0,0,1,1) << "\n";

  const arma::mat err = arma::repelem(err1, 1, n_factors) %
    arma::repmat(err2, 1, n_returns);
  Rcpp::Rcout << "err = " << err.submat(0,0,2,2) << "\n";

  W = err.t() * err;

  // for (unsigned int obs = 0; obs < n_observations; ++obs) {
  //
  //   const arma::mat err = arma::kron(err1, err2.row(obs));
  //   W += err.t() * err;
  //
  // }
  Rcpp::Rcout << "W = " << W.submat(0,0,2,2) << "\n";

  // svd of theta
  arma::mat U(n_factors, n_factors);
  arma::mat V(n_returns, n_returns);
  arma::vec sv(n_factors);

  arma::svd(U, sv, V, theta);

  Rcpp::Rcout << "U = " << U << "\n";
  Rcpp::Rcout << "V = " << V << "\n";

  arma::mat output(2, n_factors);

  // test for rank q = 0, ..., n_factors - 1
  for (unsigned int q = 0; q < n_factors; q++) {

    if (q != 1) continue;
    output.col(q) = KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
      theta, U, V, W, q, n_observations
    );

  }

  // indices of p-values above `level`
  const arma::uvec idx_accept = arma::find(output.row(1) > level);

  // the estimate of the rank is the first value `q` with associated
  // p-value above `level`. If there is no p-value above `level`, the rank is
  // set to `n_factors`
  const unsigned int rank = idx_accept.n_elem > 0 ?
    arma::min(idx_accept) : n_factors;

  return Rcpp::List::create(
    Rcpp::Named("rank") = rank,
    Rcpp::Named("ranks") = arma::regspace(0, n_factors - 2),
    Rcpp::Named("statistics") = output.row(0),
    Rcpp::Named("pvalues") = output.row(1)
  );

}

/////////////////////////////////////////////////////
///// BetaRankChenFang2019StatisticAndPvalueCpp /////

arma::vec2 BetaRankChenFang2019StatisticAndPvalueCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const unsigned int n_bootstrap,
  const double level_kp_test
) {

  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;

  if (n_factors >= n_returns) Rcpp::stop("Chen Fang 2019 test: n_factors must be < n_returns");

  const unsigned int n_observations = returns.n_rows;

  const arma::mat beta = ScaledFactorLoadingsCpp(
    returns,
    factors
  ).t();

  // SVD of beta
  arma::mat U(n_returns, n_returns);
  arma::mat V(n_factors, n_factors);
  arma::vec sv(n_factors);

  arma::svd(U, sv, V, beta);

  arma::vec2 output;

  // test statistic
  output(0) = n_observations * sv(n_factors - 1) * sv(n_factors - 1);

  // If `level_kp_test <= 0`, the initial rank estimator is taken
  // to be the number of singular values above `n_observations^(-1/3)`.
  // If `level_kp_test > 0`, the initial rank estimator is based on
  // the iterative Kleibergen Paap 2006 test with `level = level_kp_test`.
  const unsigned int rank_estimate = level_kp_test > 0. ?
  IterativeKleibergenPaap2006BetaRankTestCpp(
      returns,
      factors,
      level_kp_test
    )["rank"] :
    arma::sum(sv >= std::pow(n_observations, -1./3));

  // if full rank, set p-value to 0
  if (rank_estimate == n_factors) {
    output(1) = 0.;
    return output;
  }

  const arma::mat U22 = U.tail_cols(n_returns - rank_estimate);
  const arma::mat V22 = V.tail_cols(n_factors - rank_estimate);
  arma::vec min_sv_boot(n_bootstrap);

  const double sqrt_n_obs = std::sqrt(n_observations);

  for (unsigned int boot = 0; boot < n_bootstrap; ++boot) {

    // bootstrap indices
    const arma::uvec boot_indices = arma::randi<arma::uvec>(
      n_observations,
      arma::distr_param(0, n_observations - 1)
    );

    // bootstrap beta estimate
    const arma::mat beta_boot = ScaledFactorLoadingsCpp(
      returns.rows(boot_indices),
      factors.rows(boot_indices)
    ).t();

    // bootstrap minimum singular values
    min_sv_boot(boot) = arma::min(arma::svd(
      sqrt_n_obs * U22.t() * (beta_boot - beta) * V22
    ));

  }

  // p-value
  output(1) = (double)arma::sum(min_sv_boot % min_sv_boot >= output(0)) / n_bootstrap;

  return output;

}

///////////////////////////////////
///// ScaledFactorLoadingsCpp /////

arma::mat ScaledFactorLoadingsCpp(
  const arma::mat& returns,
  const arma::mat& factors
) {

  const arma::mat fac_t_fac = factors.t() * factors;

  const arma::mat pi = arma::solve(
    fac_t_fac, factors.t() * returns,
    arma::solve_opts::likely_sympd
  );

  return arma::chol(fac_t_fac) * arma::solve(
    arma::chol(returns.t() * returns), pi,
    arma::solve_opts::likely_sympd
  ).t();

}
