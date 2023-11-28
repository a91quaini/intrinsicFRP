// Author: Alberto Quaini

#include "hac_covariance.h"

//////////////////////////////////////////
//////// HACCovarianceMatrixCpp //////////

arma::mat HACCovarianceMatrixCpp(
  arma::mat& series,
  const bool prewhite
) {

  // Compute the number of observations from the input series
  const unsigned int n_observations = series.n_rows;

  // Determine the number of lags using the Newey-West plug-in procedure
  // as long as n_observations > 5. Otherwise set n_lags to 0
  const unsigned int n_lags =  n_observations > 5 ?
    std::floor(
      4 * std::pow(static_cast<double>(0.01 * n_observations), 2.0 / 9.0)
    ) : 0;

  // If prewhite is selected, fit an AR(1) to prewhite the series
  arma::mat coefficients;
  if (prewhite) HACPrewhiteningCpp(series, coefficients);

  // Initialize the covariance matrix with the outer product of the series
  arma::mat hac_covariance = series.t() * series /
    static_cast<double>(n_observations);

  // Vector of matrices to hold the transposed lagged series for efficiency
  std::vector<arma::mat> transposed_lagged_series(n_lags);
  for (int lag = 1; lag <= n_lags; ++lag) {

    // Store the transposed lagged series to avoid recomputing
    transposed_lagged_series[lag - 1] = series.tail_rows(n_observations - lag).t();

  }

  // Iteratively update the covariance matrix for each lag
  for (int lag = 1; lag <= n_lags; ++lag) {

    // Calculate the weight for the current lag
    const double lag_weight = 1. - static_cast<double>(lag) / (n_lags + 1.0);

    // Add the weighted lagged covariance to the estimate
    hac_covariance += lag_weight * (transposed_lagged_series[lag - 1] *
      series.head_rows(n_observations - lag)) /
      static_cast<double>(n_observations);

  }

  // If prewhite was selected, revert it
  if (prewhite) HACRevertPrewhiteningCpp(coefficients, hac_covariance);

  // Scale the covariance matrix by the number of observations
  return hac_covariance;

}

////////////////////////////////////////
//////// HACStandardErrorsCpp //////////

arma::vec HACStandardErrorsCpp(
  arma::mat& series,
  const bool prewhite
) {

  const unsigned int n_observations = series.n_rows;
  const unsigned int n_series = series.n_cols;

  // Determine the number of lags using the Newey-West plug-in procedure
  // as long as n_observations > 5. Otherwise set n_lags to 0
  const unsigned int n_lags =  n_observations > 5 ?
    std::floor(
      4 * std::pow(static_cast<double>(0.01 * n_observations), 2.0 / 9.0)
    ) : 0;

  // If prewhite is selected, fit an AR(1) to prewhite each column of the series
  arma::vec coefficients(n_series);
  if (prewhite) HACPrewhiteningCpp(series, coefficients);

  // Compute the diagonal elements of the covariance matrix for lag 0
  arma::rowvec diagonal_of_covariance = arma::sum(arma::square(series), 0) /
    static_cast<double>(n_observations);

  // Loop to compute the standard errors for lags >= 1
  for (int lag = 1; lag <= n_lags; ++lag) {

    // set the lag weight
    const double lag_weight = 1. - lag / (n_lags + 1.);

    // Loop over observations and update the diagonal with lagged contributions
    for (int obs = lag; obs < n_observations; ++obs) {

      diagonal_of_covariance += 2. * lag_weight * (
        series.row(obs - lag) % series.row(obs)
      ) / static_cast<double>(n_observations);

    }

  }

  // If prewhite was selected, revert it
  if (prewhite) HACRevertPrewhiteningCpp(coefficients, diagonal_of_covariance);

  // Return the standard errors as the square root of the marginal variances
  return arma::sqrt(diagonal_of_covariance).t();

}

////////////////////////////////////////
//////// HACStandardErrorsCpp //////////

double HACVarianceCpp(
  arma::vec& series,
  const bool prewhite
) {

  // Compute the number of observations in the series
  const unsigned int n_observations = series.n_elem;

  // Determine the number of lags using the Newey-West plug-in procedure
  // as long as n_observations > 5. Otherwise set n_lags to 0
  const unsigned int n_lags =  n_observations > 5 ?
    std::floor(
      4 * std::pow(static_cast<double>(0.01 * n_observations), 2.0 / 9.0)
    ) : 0;

  // If prewhite was selected, fit an AR(1) to prewhite the series
  double coefficient = 0.0;
  if (prewhite) HACPrewhiteningCpp(series, coefficient);

  // Initialize variance with the sum of squared observations (lag 0)
  double variance = arma::dot(series, series) /
    static_cast<double>(n_observations);

  // Loop over each lag to update the variance estimate
  for (int lag = 1; lag <= n_lags; ++lag) {
    // Calculate the Newey-West weight for the current lag
    const double lag_weight = 1. - static_cast<double>(lag) / (n_lags + 1.0);

    // Sum the weighted cross products of the series with its lagged self
    variance += 2 * lag_weight * arma::dot(
      series.tail(n_observations - lag),
      series.head(n_observations - lag)
    ) / static_cast<double>(n_observations);
  }

  // Adjust the variance using the prewhitening coefficient if prewhite is true
  if (prewhite) HACRevertPrewhiteningCpp(coefficient, variance);

  // Return the variance.
  return variance;

}

/////////////////////////////////////////////
///////// Prewhitening /////////////////////

// Matrix series
void HACPrewhiteningCpp(arma::mat& series, arma::mat& coefficients) {

  // Store the number of observations
  const unsigned int n_observations = series.n_rows;

  // Store the first n_observations - 1 lags
  const arma::mat head = series.head_rows(n_observations - 1);

  // Store the first n_observations - 1 lags
  const arma::mat tail = series.tail_rows(n_observations - 1);

  // compute the coefficients of an AR(1)
  coefficients = arma::solve(head, tail);

  // Prewhite the series using the coefficients
  series = arma::join_vert(
    series.row(0),
    tail - head * coefficients
  );

}

// Matrix series (marginal pre-whitening)
void HACPrewhiteningCpp(arma::mat& series, arma::vec& coefficients) {

  // Store the number of observations
  const unsigned int n_observations = series.n_rows;

  // Loop over the colummns of the series
  for (int col = 0; col < series.n_cols; ++col) {

    // Store the first n_observations - 1 lags
    const arma::vec head = series.col(col).head(n_observations - 1);

    // Store the last n_observations - 1 lags
    const arma::vec tail = series.col(col).tail(n_observations - 1);

    // compute the coefficients of an AR(1)
    coefficients(col) = arma::dot(head, tail) / arma::dot(head, head);

    const arma::vec first_obs = series(arma::span(0), 0);

    // Prewhite the series using the coefficients
    series.col(col) = arma::join_vert(
      series(arma::span(0), col),
      tail - head * coefficients(col)
    );

  }

}

// Scalar series
void HACPrewhiteningCpp(arma::vec& series, double coefficient) {

  // Store the number of observations
  const unsigned int n_observations = series.n_rows;

  // Store the first n_observations - 1 lags
  const arma::vec head = series.head(n_observations - 1);

  // Store the last n_observations - 1 lags
  const arma::vec tail = series.tail(n_observations - 1);

  // Compute the coefficients of an AR(1)
  coefficient = arma::dot(head, tail) / arma::dot(head, head);

  // Prewhite the series using the coefficients
  series = arma::join_vert(
    series(arma::span(0)),
    tail - head * coefficient
  );

}

// Matrix series
void HACRevertPrewhiteningCpp(
  const arma::mat& coefficients,
  arma::mat& hac_covariance
) {

  // Store the number of variables
  const unsigned int n_variables = coefficients.n_rows;

  // Compute the inverse of the identity -  coefficients.t()
  arma::mat inv_temp;
  try {
    // Try arma::inv_sympd
    inv_temp = arma::inv(arma::eye(n_variables, n_variables) - coefficients.t(), arma::inv_opts::no_ugly);

  } catch (const std::runtime_error&) {
    // Fallback to generalized inverse
    inv_temp = pinv(
      arma::eye(n_variables, n_variables) - coefficients.t()
    );

  }

  // Adjust the HAC covariance matrix using the inverse of the pre-whitening transformation
  hac_covariance = inv_temp * hac_covariance * inv_temp.t();

}

// Matrix series (marginal pre-whitening)
void HACRevertPrewhiteningCpp(
  const arma::vec& coefficients,
  arma::rowvec& hac_covariance
) {

  // Loop over the elements of hac_covariance
  for (int idx = 0; idx < hac_covariance.n_elem; ++idx) {

    // Compute 1 - the coefficient(idx)
    const double temp = 1. - coefficients(idx);

    // Adjust the HAC covariance element using the inverse of the prewhitening transformation
    hac_covariance(idx) /= (temp * temp);

  }

}

// Scalar series
void HACRevertPrewhiteningCpp(
  const double coefficient,
  double hac_covariance
) {

  // Compute 1 - the coefficient
  const double temp = 1. - coefficient;

  // Adjust the HAC covariance using the inverse of the prewhitening transformation
  hac_covariance /= (temp * temp);

}
