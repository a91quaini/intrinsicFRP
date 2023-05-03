// Author: Alberto Quaini

#include "hac_standard_errors.h"

//////////////////////////////////
///// HACStandardErrorsCpp ///////

arma::vec HACStandardErrorsCpp(const arma::mat& series) {

  const unsigned int n_observations = series.n_rows;

  const unsigned int n_lags = 4 * std::pow(n_observations/100, 2./9.);

  // lag 0
  arma::rowvec diagonal_of_covariance = arma::diagvec(series.t() * series).t();

  // lags >= 1
  for (unsigned int lag = 1; lag <= n_lags; ++lag) {

    const double lag_weight = 1. - lag / (n_lags + 1.);

    for (unsigned int obs = lag; obs < n_observations; ++obs) {

      diagonal_of_covariance += 2. * lag_weight * (
        series.row(obs - lag) % series.row(obs)
      );

    }

  }

  return arma::sqrt(diagonal_of_covariance).t() / n_observations;

}
