// Author: Alberto Quaini

#include "identification_tests.h"
#include "utils.h"

///////////////////////////////////////
///// ChenFang2019BetaRankTestCpp /////

Rcpp::List ChenFang2019BetaRankTestCpp(
    const arma::mat& returns, // Asset returns matrix
    const arma::mat& factors, // Factor returns matrix
    const unsigned int n_bootstrap, // Number of bootstrap samples for p-value calculation
    const double target_level_kp2006_rank_test // Target level for the KP2006 rank test
) {

  // Ensure that the number of factors is less than the number of returns
  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;
  if (n_factors >= n_returns) Rcpp::stop("Chen Fang 2019 test: n_factors must be < n_returns");

  // Number of observations is the number of rows in the returns matrix
  const unsigned int n_observations = returns.n_rows;

  // Calculate scaled factor loadings and transpose the result for later operations
  const arma::mat beta = ScaledFactorLoadingsCpp(returns, factors).t();

  // Perform Singular Value Decomposition (SVD) on the transposed beta matrix
  arma::mat U(n_returns, n_returns);
  arma::mat V(n_factors, n_factors);
  arma::vec sv(n_factors);
  arma::svd(U, sv, V, beta);

  // Calculate the test statistic using the smallest singular value
  const double statistic = static_cast<double>(n_observations) *
    std::pow(sv(n_factors - 1), 2);

  // Estimate the rank of beta matrix based on the target level for KP2006 test
  // Choose method based on the target level specified
  const unsigned int rank_estimate = target_level_kp2006_rank_test > 0. ?
    IterativeKleibergenPaap2006BetaRankTestCpp(
      returns, factors, target_level_kp2006_rank_test
    )["rank"] :
    arma::sum(sv >= std::pow(static_cast<double>(n_observations), -1./4));

  // If estimated rank is full, set the p-value to 0 as the test statistic will not follow the desired distribution
  if (rank_estimate == n_factors) {

    return Rcpp::List::create(
      Rcpp::Named("statistic") = statistic,
      Rcpp::Named("p-value") = 0.
    );

  }

  // Extract the bottom right blocks of U and V matrices for bootstrap calculation
  const arma::mat U22 = U.tail_cols(n_returns - rank_estimate);
  const arma::mat V22 = V.tail_cols(n_factors - rank_estimate);

  // Vector to store the minimum singular values from bootstrap samples
  arma::vec min_sv_boot(n_bootstrap);

  // Square root of the number of observations, used in bootstrapping
  const double sqrt_n_obs = std::sqrt(static_cast<double>(n_observations));

  // Bootstrap loop to calculate p-value
  for (unsigned int boot = 0; boot < n_bootstrap; ++boot) {

    // Generate random bootstrap indices
    const arma::uvec boot_indices = arma::randi<arma::uvec>(
      n_observations,
      arma::distr_param(0, n_observations - 1)
    );

    // Calculate bootstrap beta estimate
    const arma::mat beta_boot = ScaledFactorLoadingsCpp(
      returns.rows(boot_indices),
      factors.rows(boot_indices)
    ).t();

    // Calculate minimum singular value for the bootstrap sample
    min_sv_boot(boot) = arma::min(arma::svd(
      sqrt_n_obs * U22.t() * (beta_boot - beta) * V22
    ));

  }

  // Calculate the p-value as the proportion of bootstrap samples where
  // the squared minimum singular value is greater than or equal to the test statistic
  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("p-value") = (double)arma::sum(
      min_sv_boot % min_sv_boot >= statistic
    ) / static_cast<double>(n_bootstrap)
  );

}


///////////////////////////////////////////////////////////////
///// KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp /////

arma::vec2 KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
  const arma::mat& theta_vectorised,
  const arma::mat& U,
  const arma::mat& V,
  const arma::mat& W,
  const unsigned int q
) {

  // Determine the number of factors and returns from the dimensions of U and V
  const unsigned int n_factors = U.n_rows;
  const unsigned int n_returns = V.n_cols;

  // Extract the bottom-right submatrices after the hypothesized rank (q)
  // according to equation (3) in the KP 2006 paper
  const arma::mat U22 = U.submat(q, q, n_factors - 1, n_factors - 1);
  const arma::mat V22 = V.submat(q, q, n_returns - 1, n_returns - 1);

  ///////////////////////////////////////////
  // // construct of symmetric square-roots of matrices (U22 U22') and (V22 V22')
  // // need to clamp the eigenvalues to avoid spurious complex numbers
  // // const arma::mat sqrt_U22 = arma::sqrtmat_sympd(U22 * U22.t());
  // // const arma::mat sqrt_V22 = arma::sqrtmat_sympd(V22 * V22.t());
  // arma::vec U22_eval(n_factors - q);
  // arma::vec V22_eval(n_factors - q);
  // arma::mat U22_evec(n_factors - q, n_factors - q);
  // arma::mat V22_evec(n_factors - q, n_factors - q);
  //
  // eig_sym(U22_eval, U22_evec, U22 * U22.t());
  // eig_sym(V22_eval, V22_evec, V22 * V22.t());
  //
  // const arma::vec sqrt_U22_eval = arma::sqrt(U22_eval);
  // const arma::vec sqrt_V22_eval = arma::sqrt(V22_eval);
  //
  // const arma::mat sqrt_U22 = U22_evec * arma::diagmat(
  //   arma::clamp(sqrt_U22_eval, 0., sqrt_U22_eval.max())
  // ) * U22_evec.t();
  // const arma::mat sqrt_V22 = V22_evec * arma::diagmat(
  //   arma::clamp(sqrt_V22_eval, 0., sqrt_V22_eval.max())
  // ) * V22_evec.t();
  //////////////////////////////////////////

  // Construct square roots directly using eigen decomposition results.
  const arma::mat sqrt_U22 = arma::sqrtmat_sympd(U22 * U22.t());
  const arma::mat sqrt_V22 = arma::sqrtmat_sympd(V22 * V22.t());

  // Compute the A_qperp and B_qperp matrices according to equation (12) in the KP 2006 paper
  const arma::mat A_qperp = U.tail_cols(n_factors - q) * SolveSympd(
    U22,
    sqrt_U22
  );
  const arma::mat B_qperp = sqrt_V22 * SolveSympd(
    V22.t(),
    V.tail_cols(n_returns - q).t()
  );

  // Compute the Kronecker product of B_qperp and A_qperp transposed
  const arma::mat kron_BA_qperp = arma::kron(B_qperp, A_qperp.t());

  // Construct the lambda_q vector, which represents the minimized eigenvalues
  // of the concentrated parameter space (below equation (21) in KP 2006)
  const arma::vec lambda_q = kron_BA_qperp * theta_vectorised;

  // Compute Omega_q, the covariance matrix of lambda_q
  const arma::mat Omega_q = kron_BA_qperp * W * kron_BA_qperp.t();

  // Calculate the test statistic, which does not require multiplication by n_observations
  // This is based on the quadratic form of lambda_q and Omega_q
  const double statistic = arma::dot(lambda_q, SolveSympd(
    Omega_q,
    lambda_q
  ));

  // Calculate the p-value using the Chi-squared distribution with degrees of freedom
  // given by the product of the excess dimensions of factors and returns over hypothesized rank q
  const double pvalue = R::pchisq(
    statistic, (n_factors - q) * (n_returns - q),
    false,
    false
  );

  // Return a 2-element vector containing the test statistic and its p-value
  return arma::vec2{statistic, pvalue};

}

//////////////////////////////////////////////////////
///// IterativeKleibergenPaap2006BetaRankTestCpp /////

Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const double target_level
) {

  // Obtain the number of returns and factors from the size of the matrices
  const unsigned int n_returns = returns.n_cols;
  const unsigned int n_factors = factors.n_cols;

  // Ensure the number of factors is less than the number of returns
  if (n_factors >= n_returns) Rcpp::stop("Kleibergen Paap test: n_factors must be < n_returns");

  // Perform Cholesky decomposition of factors.t() * factors.
  // This is the transpose of matrix G in equation (16) in Kleibergen-Paap (2006)
  const arma::mat L_fac = arma::chol(factors.t() * factors, "lower");

  // Compute the Cholesky decomposition of returns.t() * returns.
  // This is the transpose of matrix F_t_inv in equation (16) in Kleibergen-Paap (2006)
  const arma::mat U_ret = arma::chol(returns.t() * returns);

  // Solve for the coefficients matrix (pi) using the factors and returns matrices
  // This follows equation (16) in Kleibergen-Paap (2006)
  const arma::mat pi = arma::solve(arma::trimatl(L_fac), factors.t() * returns);

  // Compute matrix theta of factor loadings that is proportional to matrix of
  // t-statistics of the least squares estimator.
  // That is, it is invariant to invertible transformations of the data that are
  // identical over all observations.
  const arma::mat theta = pi * arma::inv(arma::trimatu(U_ret));

  // Compute centred returns
  const arma::mat returns_centred = returns.each_row() - arma::mean(returns);

  // Calculate the residuals and then apply scaling by F_t_inv and G
  const arma::mat residuals = returns_centred - factors * pi;
  const arma::mat err1 = arma::solve(arma::trimatl(U_ret.t()), residuals.t()).t();
  const arma::mat err2 = arma::solve(arma::trimatl(L_fac), factors.t()).t();

  // Element-wise multiplication of errors, repeated across the correct dimensions to match W
  const arma::mat err = arma::repelem(err1, 1, n_factors) % arma::repmat(err2, 1, n_returns);

  // Compute the White covariance matrix with the outer product of errors
  const arma::mat W = err.t() * err;

  // Perform singular value decomposition (SVD) on theta
  arma::mat U(n_factors, n_factors);
  arma::mat V(n_returns, n_returns);
  arma::vec sv(n_factors);
  arma::svd(U, sv, V, theta);

  arma::mat output(2, n_factors); // Matrix to hold the test statistics and p-values

  // Iteratively test for each possible rank q from 0 to n_factors - 1
  for (unsigned int q = 0; q < n_factors; q++) {

    // Compute the test statistic and p-value for the current rank q
    output.col(q) = KleibergenPaap2006BetaRankTestStatisticAndPvalueCpp(
      arma::vectorise(theta), U, V, W, q
    );

  }

  // Find the indices of p-values that are above the significance level, adjusted for multiple testing
  const arma::uvec idx_accept = arma::find(output.row(1) > target_level / factors.n_cols);

  // Estimate the rank as the smallest q with a p-value above the significance level
  // If no p-value is above the level, the estimated rank is set to n_factors
  const unsigned int rank = idx_accept.empty() ? n_factors : arma::min(idx_accept);

  // Return a list with the estimated rank, the range of q tested, and the test statistics and p-values
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

  // Perform Cholesky decomposition of factors.t() * factors.
  const arma::mat L_fac = arma::chol(factors.t() * factors, "lower");

  // Compute the Cholesky decomposition of returns.t() * returns.
  const arma::mat U_ret = arma::chol(returns.t() * returns);

  // Compute pi as the solution x to L_fac * x = factors.t() * returns.
  const arma::mat pi = arma::solve(arma::trimatl(L_fac), factors.t() * returns);

  // Compute matrix of factor loadings that is proportional to matrix of
  // t-statistics of the least squares estimator.
  // That is, it is invariant to invertible transformations of the data that are
  // identical over all observations.
  return pi * arma::inv(arma::trimatu(U_ret));

}
