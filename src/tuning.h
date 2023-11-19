// Author: Alberto Quaini

#ifndef TUNING_H
#define TUNING_H

#include <RcppArmadillo.h>

// Function for internal use
// Computes the Generalized Cross Validation (GCV) score for each penalty parameter value for Oracle Tradable Factor Risk Premia (TFRP).
// It evaluates the fit of the model to the data, balancing the trade-off between model complexity and accuracy.
// The function considers optional arguments for scaling, model identification checks, and a target level for the Kleibergen Paap 2006 rank test.
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
// Computes the Cross Validation (CV) score for Oracle Tradable Factor Risk Premia (TFRP) for each penalty parameter value.
// It helps to select the optimal penalty parameter by evaluating the model's predictive performance on separate folds of the data set.
arma::vec CVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type = 'c',
  const unsigned int n_folds = 5
);

// Function for internal use
// Computes the Rolling Validation (RV) score for Oracle Tradable Factor Risk Premia (TFRP) for each penalty parameter value.
// It helps to select the optimal penalty parameter by evaluating the model's predictive performance  over rolling training and testing windows
// thereby improving the model's stability and predictive ability over time.
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
// Computes the prediction error of the model on the test data.
// This function calculates the discrepancy between the predicted and actual returns,
// which helps in assessing the model's accuracy.
double ComputePredictionErrorCpp(
  const arma::vec& oracle_frp_train,
  const arma::mat& covariance_factor_returns_train,
  const arma::mat& variance_returns_train,
  const arma::vec& mean_returns_test,
  const arma::vec& variance_returns_test
);

// Function for internal use
// Computes the Generalized Cross Validation (GCV) score. This function is utilized for model selection,
// providing a measure of fit that penalizes model complexity to prevent overfitting.
double ComputeGCVScoreCpp(
  const arma::vec& oracle_tfrp,
  const arma::uvec& idx_selected,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double scaling
);

#endif
