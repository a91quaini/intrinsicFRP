// Author: Alberto Quaini

#ifndef HAC_STANDARD_ERRORS_H
#define HAC_STANDARD_ERRORS_H

#include <RcppArmadillo.h>

// For internal use
// Computes the Heteroskedasticity and Autocorrelation robust standard errors
// of n_variables random variables based on a
// (n_observations x n_variables)-dimensional matrix of realizations.
// It uses the Newey-West (1994) plug-in procedure to select the number of
// relevant lags, i.e., n_lags = 4 * (n_observations/100)^(2/9).
arma::vec HACStandardErrorsCpp(const arma::mat& series);

#endif
