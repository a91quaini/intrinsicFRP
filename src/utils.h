// Author: Alberto Quaini

#ifndef UTILS_H
#define UTILS_H

#include <RcppArmadillo.h>

// Set template T, representing both arma::mat or arma::vec in subsequent functions.
template<typename T>

// Function for internal use.
// Solve a system of linear equations
// A * solution = b
// where A is likely symmetric and positive definite,
// and b is a vector or a matrix.
arma::mat SolveSympd(const arma::mat& A, const T& b) {

  try {
    // Try arma::solve
    return arma::solve(
      A, b, arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
    );

  } catch (const std::runtime_error&) {
    // Fallback to generalized inverse
    return arma::solve(
      A, b, arma::solve_opts::force_approx
    );

  }

}

// Function for internal use.
// Invert a symmetric and positive definite matrix A.
arma::mat InvSympd(const arma::mat& A);

#endif
