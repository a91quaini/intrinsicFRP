// Author: Alberto Quaini

#include "adaptive_weights.h"

///////////////////////////////
///// AdaptiveWeightsCpp //////

arma::vec AdaptiveWeightsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const char type
) {

  switch (type) {

    case 'c' : return AdaptiveWeightsFromMatrixCpp(arma::cor(factors, returns));

    case 'a' : return AdaptiveWeightsFromVectorCpp(
        arma::cov(factors, returns) * arma::solve(
        arma::cov(returns),
        arma::mean(returns).t(),
        arma::solve_opts::likely_sympd
      ));

    case 'b' : return AdaptiveWeightsFromMatrixCpp(arma::solve(
        arma::cov(factors),
        arma::cov(factors, returns)
      ));

    default : return arma::ones(factors.n_cols);

  }

}

arma::vec AdaptiveWeightsFromMatrixCpp(const arma::mat& matrix) {

  return 1. / arma::sum(matrix % matrix, 1);

}

arma::vec AdaptiveWeightsFromVectorCpp(const arma::vec& vector) {

  return 1. / vector % vector;

}
