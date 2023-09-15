// Author: Alberto Quaini

#ifndef TUNING_H
#define TUNING_H

#include <RcppArmadillo.h>

// Function for internal use
// Computes the Generalized Cross Validation score of the oracle tradable
// factor risk premia for each penalty parameter value.
arma::vec GCVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& oracle_frp,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const unsigned int n_observations,
  const bool gcv_scaling_n_assets = false,
  const bool gcv_identification_check = false,
  const double target_level_kp2006_rank_test = 0.05
);

// Function for internal use
// Computes the Cross Validation score of the oracle tradable factor risk
// premia for each penalty parameter value.
arma::vec CVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const unsigned int n_folds = 5
);

// Function for internal use
// Computes the Rolling Validation score of the oracle tradable factor risk
// premia for each penalty parameter value.
arma::vec RVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const unsigned int n_train_observations = 120,
  const unsigned int n_test_observations = 12,
  const unsigned int roll_shift = 12
);

// Function for internal use
// Computes the prediction error given information on the train and test data.
double ComputePredictionErrorCpp(
  const arma::vec& oracle_frp_train,
  const arma::mat& covariance_factor_returns_train,
  const arma::mat& variance_returns_train,
  const arma::vec& mean_returns_test,
  const arma::vec& variance_returns_test
);

// Function for internal use
// Computes the GCV score.
double ComputeGCVScoreCpp(
  const arma::vec& oracle_tfrp,
  const arma::uvec& idx_selected,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double scaling
);

#endif
