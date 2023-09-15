// Author: Alberto Quaini

#include "identification_tests.h"

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

  // test statistic
  const double statistic = n_observations * sv(n_factors - 1) * sv(n_factors - 1);

  // If `target_level_kp2006_rank_test <= 0`, the initial rank estimator is taken
  // to be the number of singular values above `n_observations^(-1/4)`.
  // If `target_level_kp2006_rank_test > 0`, the initial rank estimator is based on
  // the iterative Kleibergen Paap 2006 test with
  // `level = target_level_kp2006_rank_test / n_factors`.
  const unsigned int rank_estimate = target_level_kp2006_rank_test > 0. ?
  IterativeKleibergenPaap2006BetaRankTestCpp(
    returns,
    factors,
    target_level_kp2006_rank_test
  )["rank"] :
    arma::sum(sv >= std::pow(n_observations, -1./4));

  // if full rank, set p-value to 0
  if (rank_estimate == n_factors) {

    return Rcpp::List::create(
      Rcpp::Named("statistic") = statistic,
      Rcpp::Named("p-value") = 0.
    );

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

  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("p-value") = (double)arma::sum(
      min_sv_boot % min_sv_boot >= statistic
    ) / n_bootstrap
  );

}


///////////////////////////////////////////////////////////////
///// KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp /////

arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::mat& theta_vectorised,
  const arma::mat& U,
  const arma::mat& V,
  const arma::mat& W,
  const unsigned int q,
  const unsigned int n_observations
) {

  const unsigned int n_factors = U.n_rows;
  const unsigned int n_returns = V.n_cols;

  // see eq. (3) in Kleibergen Paap 2006
  const arma::mat U22 = U.submat(q, q, n_factors - 1, n_factors - 1);
  const arma::mat V22 = V.submat(q, q, n_returns - 1, n_returns - 1);

  // construct of symmetric square-roots of matrices (U22 U22') and (V22 V22')
  // need to clamp the eigenvalues to avoid spurious complex numbers
  // const arma::mat sqrt_U22 = arma::sqrtmat_sympd(U22 * U22.t());
  // const arma::mat sqrt_V22 = arma::sqrtmat_sympd(V22 * V22.t());
  arma::vec U22_eval(n_factors - q);
  arma::vec V22_eval(n_factors - q);
  arma::mat U22_evec(n_factors - q, n_factors - q);
  arma::mat V22_evec(n_factors - q, n_factors - q);

  eig_sym(U22_eval, U22_evec, U22 * U22.t());
  eig_sym(V22_eval, V22_evec, V22 * V22.t());

  const arma::vec sqrt_U22_eval = arma::sqrt(U22_eval);
  const arma::vec sqrt_V22_eval = arma::sqrt(V22_eval);

  const arma::mat sqrt_U22 = U22_evec * arma::diagmat(
    arma::clamp(sqrt_U22_eval, 0., sqrt_U22_eval.max())
  ) * U22_evec.t();
  const arma::mat sqrt_V22 = V22_evec * arma::diagmat(
    arma::clamp(sqrt_V22_eval, 0., sqrt_V22_eval.max())
  ) * V22_evec.t();

  // see eq. (12) in Kleibergen Paap 2006
  const arma::mat A_qperp = U.tail_cols(n_factors - q) * arma::solve(
    U22, sqrt_U22,
    arma::solve_opts::likely_sympd
  );
  const arma::mat B_qperp = sqrt_V22 * arma::solve(
    V22.t(), V.tail_cols(n_returns - q).t(),
    arma::solve_opts::likely_sympd
  );

  const arma::mat kron_BA_qperp = arma::kron(B_qperp, A_qperp.t());

  // see below eq. (21) in Kleibergen Paap 2006
  // we use that vec(A'lB') = (B kron A') vec(l)
  const arma::vec lambda_q = kron_BA_qperp * theta_vectorised;
  const arma::mat Omega_q = kron_BA_qperp * W * kron_BA_qperp.t();

  // test statistic: no need to multiply by n_observations
  const double statistic = arma::dot(lambda_q, arma::solve(
    Omega_q, lambda_q,
    arma::solve_opts::likely_sympd
  ));

  // p-value: probability of observing a larger value than the test statistic
  // under a Chi-square((n_factors - q) * (n_returns - q))
  const arma::vec2 output = {
    statistic,
    R::pchisq(
      statistic,
      (n_factors - q) * (n_returns - q),
      false,
      false
    )
  };

  return output;

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

  if (n_factors >= n_returns) Rcpp::stop("Kleibergen Paap test: n_factors must be < n_returns");

  const unsigned int n_observations = factors.n_rows;

  const arma::mat fac_t_fac = factors.t() * factors;

  // see below eq. (16) in Kleibergen Paap 2006
  const arma::mat pi = arma::solve(
    fac_t_fac, factors.t() * returns,
    arma::solve_opts::likely_sympd
  );

  // see eq. (16) in Kleibergen Paap 2006
  const arma::mat F_t_inv = arma::chol(returns.t() * returns);
  const arma::mat G = arma::chol(fac_t_fac);
  const arma::mat theta = G * arma::solve(
    arma::trimatl(F_t_inv.t()), pi.t()
  ).t();

  // White covariance matrix estimator with scaling matrices incorporated
  arma::mat W(n_returns * n_factors, n_returns * n_factors, arma::fill::eye);

  const arma::mat residuals = returns - factors * pi - arma::repmat(
    arma::mean(returns), n_observations, 1
  );
  const arma::mat err1 = arma::solve(
    arma::trimatl(F_t_inv.t()), residuals.t()
  ).t();
  const arma::mat err2 = arma::solve(
    arma::trimatl(G.t()), factors.t()
  ).t();

  const arma::mat err = arma::repelem(err1, 1, n_factors) %
    arma::repmat(err2, 1, n_returns);

  W = err.t() * err;

  // svd of theta
  arma::mat U(n_factors, n_factors);
  arma::mat V(n_returns, n_returns);
  arma::vec sv(n_factors);

  arma::svd(U, sv, V, theta);

  arma::mat output(2, n_factors);

  // test for rank q = 0, ..., n_factors - 1
  for (unsigned int q = 0; q < n_factors; q++) {

    // if (q != 1) continue;
    output.col(q) = KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
      arma::vectorise(theta), U, V, W, q, n_observations
    );

  }

  // indices of p-values above `level`
  const arma::uvec idx_accept = arma::find(
    output.row(1) > target_level / factors.n_cols
  );

  // the estimate of the rank is the first value `q` with associated
  // p-value above `level`. If there is no p-value above `level`, the rank is
  // set to `n_factors`
  const unsigned int rank = idx_accept.n_elem > 0 ?
    arma::min(idx_accept) : n_factors;

  return Rcpp::List::create(
    Rcpp::Named("rank") = rank,
    Rcpp::Named("q") = arma::regspace(0, n_factors - 1).t(),
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

  const arma::mat fac_t_fac = factors.t() * factors;

  const arma::mat pi = arma::solve(
    fac_t_fac, factors.t() * returns,
    arma::solve_opts::likely_sympd
  );

  return arma::chol(fac_t_fac) * arma::solve(
      arma::trimatl(arma::chol(returns.t() * returns).t()),
      pi.t()
  ).t();

}
