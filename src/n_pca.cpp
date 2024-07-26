// Author: Alberto Quaini

#include "n_pca.h"
#include "utils.h"

////////////////////////////////////////////
/////// NPCA_GiglioXiu2021Cpp //////////////

unsigned int NPCA_GiglioXiu2021Cpp(
  const arma::vec& e_vals,
  const unsigned int n_assets,
  const unsigned int n_observations,
  unsigned int p_max
) {

  // minimum between n_assets and n_observations
  const unsigned int min_NT = std::min(n_assets, n_observations);

  // adjust p_max if it is too large
  if (p_max > min_NT) {
    p_max = min_NT - 1;
  }

  // compute the penalty term
  const double scaling = 0.5 * arma::median(e_vals);
  const double penalty = scaling * (std::log(n_assets) + std::log(n_observations)) *
    (std::pow(n_assets, -0.5) + std::pow(n_observations, -0.5));

  // compute the criterion for each index
  arma::vec criterion(p_max);
  for (unsigned int i = 0; i < p_max; ++i) {
    criterion(i) = e_vals(i) + (i + 1) * penalty;
  }

  // return the index that minimizes the criterion
  return criterion.index_min();

}

////////////////////////////////////////////////
/////// NPCA_AhnHorenstein2013Cpp //////////////

Rcpp::List NPCA_AhnHorenstein2013Cpp(
  const arma::vec& evals,
  unsigned int n_max
) {

  // set the number of eigenvalues
  const unsigned int n_evals = evals.n_elem;

  // if n_max <= zero or >= n_evals, set it to n_evals - 1
  if (n_max <= 0 || n_max >= n_evals) {
    n_max = n_evals - 1;
  }

  const arma::uvec top = arma::regspace<arma::uvec>(0, n_max - 1);
  const arma::uvec bottom = arma::regspace<arma::uvec>(1, n_max);

  // calculate the eigenvalue ratios and the index where the ratio is maximized
  const arma::vec evals_ratio = evals(top) / evals(bottom);
  const unsigned int er = evals_ratio.index_max() + 1; // +1 to match R's 1-based index

  // calculate the sum of the eigenvalues after each index
  arma::vec evals_revcumsum = arma::reverse(arma::cumsum(arma::reverse(evals)));
  evals_revcumsum.shed_row(n_evals - 1); // remove the last element
  evals_revcumsum.insert_rows(n_evals - 1, arma::zeros(1)); // add zero at the end

  // divide the eigenvalues by the sum
  const arma::vec evals_tr = evals / evals_revcumsum;
  const arma::vec growth_ratio = arma::log(
    1 + evals_tr(top)) / arma::log(1 + evals_tr(bottom)
  );
  const unsigned int gr = growth_ratio.index_max() + 1; // +1 to match R's 1-based index

  return Rcpp::List::create(
    Rcpp::Named("er") = er,
    Rcpp::Named("gr") = gr
  );
}
