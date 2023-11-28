// Author: Alberto Quaini

#include "utils.h"

///////////////////////////////
/////// InvSympd //////////////

arma::mat InvSympd(const arma::mat& A) {

  try {
    // Try arma::inv_sympd
    return arma::inv_sympd(A, arma::inv_opts::no_ugly);

  } catch (const std::runtime_error&) {
    // Fallback to generalized inverse
    return  arma::solve(
      A, arma::eye(A.n_rows, A.n_rows), arma::solve_opts::force_approx
    );

  }

}
