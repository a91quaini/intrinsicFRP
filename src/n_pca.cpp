// Author: Alberto Quaini

#include "n_pca.h"
#include "utils.h"

////////////////////////////////////////////
/////// NPCA_GiglioXiu2021Cpp //////////////

unsigned int NPCA_GiglioXiu2021Cpp(
  const arma::vec& evals,
  const unsigned int n_assets,
  const unsigned int n_observations,
  unsigned int n_max
) {

  const unsigned int min_NT = std::min(n_assets, n_observations);

  // if n_max <= zero or >= n_assets, set it to n_assets - 1
  if (n_max <= 0 || n_max >= min_NT) {
    n_max = min_NT - 1;
  }

  const arma::vec e_vals = evals.head(n_max);

  // compute the penalty term
  const double scaling = 0.5 * arma::median(e_vals);
  const double penalty = scaling * (std::log(n_assets) + std::log(n_observations)) *
    (std::pow(n_assets, -0.5) + std::pow(n_observations, -0.5));

  // compute the criterion for each index
  const arma::vec criterion = e_vals + penalty * arma::regspace(1, n_max);

  // return the number of pca
  return criterion.index_min() + 1;

}

////////////////////////////////////////////////
/////// NPCA_AhnHorenstein2013Cpp //////////////

Rcpp::List NPCA_AhnHorenstein2013Cpp(
  const arma::vec& evals,
  const unsigned int n_assets,
  const unsigned int n_observations,
  unsigned int n_max
) {

  const unsigned int min_NT = std::min(n_assets, n_observations);

  // if n_max <= zero or >= n_assets, set it to n_assets - 1
  if (n_max <= 0 || n_max >= min_NT) {
    n_max = min_NT - 1;
  }

  const arma::vec e_vals = evals.head(n_max);
  const arma::uvec top = arma::regspace<arma::uvec>(0, n_max - 2);
  const arma::uvec bottom = arma::regspace<arma::uvec>(1, n_max - 1);

  // calculate the eigenvalue ratios and the index where the ratio is maximized
  const arma::vec evals_ratio = e_vals(top) / e_vals(bottom);
  const unsigned int er = evals_ratio.index_max() + 1; // n_pca

  // calculate the sum of the eigenvalues after each index
  arma::vec evals_revcumsum = arma::reverse(arma::cumsum(arma::reverse(e_vals)));
  evals_revcumsum.shed_row(e_vals.n_elem - 1); // remove the last element
  evals_revcumsum.insert_rows(e_vals.n_elem - 1, arma::zeros(1)); // add zero at the end

  // divide the eigenvalues by the sum
  const arma::vec evals_tr = e_vals / evals_revcumsum;
  const arma::vec growth_ratio = arma::log(
    1 + evals_tr(top)) / arma::log(1 + evals_tr(bottom)
  );
  const unsigned int gr = growth_ratio.index_max() + 1; // n_pca

  return Rcpp::List::create(
    Rcpp::Named("er") = er,
    Rcpp::Named("gr") = gr
  );
}
