// Author: Alberto Quaini

#include "adaptive_weights.h"
#include "utils.h"

///////////////////////////////
///// AdaptiveWeightsCpp //////

arma::vec AdaptiveWeightsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const char type
) {

  // Choose the type of adaptive weights based on input 'type'
  switch (type) {

  case 'c' :
    // Type 'c': Compute weights based on the correlation matrix between factors and returns
    return AdaptiveWeightsFromMatrixCpp(arma::cor(factors, returns));

  case 'a' :
    // Type 'a': Compute weights based on the first-step intrinsic factor risk premia estimates
    return AdaptiveWeightsFromVectorCpp(
      arma::cov(factors, returns) * SolveSympd(
        arma::cov(returns),
        arma::mean(returns).t()
    ));

  case 'b' :
    // Type 'b': Compute weights based on the matrix of factors regression coefficients on test asset excess returns
    return AdaptiveWeightsFromMatrixCpp(SolveSympd(
      arma::cov(factors),
      arma::cov(factors, returns)
    ));

  default :
    // Other types: Return unit vector (equal weights)
    return arma::ones(factors.n_cols);

  }

}

arma::vec AdaptiveWeightsFromMatrixCpp(const arma::mat& matrix) {

  // Return the inverse of the sum of squares of each row in the matrix
  return 1. / arma::sum(arma::square(matrix), 1);

}

arma::vec AdaptiveWeightsFromVectorCpp(const arma::vec& vector) {

  // Return the inverse of the element-wise square of the vector
  return 1. / arma::square(vector);

}
