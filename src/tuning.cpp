// Author: Alberto Quaini

#include "tuning.h"
#include "tfrp.h"
#include "oracle_tfrp.h"
#include "adaptive_weights.h"
#include "identification_tests.h"
#include "utils.h"

/////////////////////////////////
///// GCVTuningOracleTFRPCpp /////

arma::vec GCVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& oracle_tfrp,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const unsigned int n_observations,
  const bool gcv_scaling_n_assets,
  const bool gcv_identification_check,
  const double target_level_kp2006_rank_test
) {

  // Initialize the score vector with values equal to the inner product of the mean of returns (no model case)
  arma::vec score(oracle_tfrp.n_cols, arma::fill::value(
      arma::accu(arma::square(mean_returns))
  ));

  // Determine the scaling factor based on the option gcv_scaling_n_assets
  const double scaling = gcv_scaling_n_assets ?
  std::sqrt(static_cast<double>(mean_returns.n_elem)) / static_cast<double>(n_observations) :
    1.0 / static_cast<double>(n_observations);

  //// Prepare variables for identification check.
  // If gcv_identification_check is selected, set is_model_identified to false,
  // otherwise set it to true, so that we never check for model identification
  // in the loop below.
  bool is_model_identified = !gcv_identification_check;

  // Store previous_model_size as n_factors + 1,
  // so that we always check the largest (first) model for identification.
  unsigned int previous_model_size = factors.n_cols + 1;
  ////

  // Iterate over each column of oracle_tfrp
  for (unsigned int par = 0; par < oracle_tfrp.n_cols; ++par) {

    // Store the indices of the nonzero elements in the appropriate oracle_tfrp
    const arma::uvec idx_selected = arma::find(oracle_tfrp.col(par));

    // Skip iteration if no factors are selected
    if (idx_selected.empty()) continue;

    // If the previous model was not identified,
    // check if the current model is identified using Kleibergen Paap 2006 Beta Rank Test.
    if (!is_model_identified) {

      // Skip iteration if the current model's size is as large
      // as the previous one (even if it's composed differently).
      if (idx_selected.n_elem >= previous_model_size) continue;

      // Update the model's size.
      previous_model_size = idx_selected.n_elem;

      // Compute the Kleibergen Paap 2006 Beta Rank Test
      const unsigned int rank_estimate = IterativeKleibergenPaap2006BetaRankTestCpp(
        returns,
        factors.cols(idx_selected),
        target_level_kp2006_rank_test
      )["rank"];

      // The model is identified if the rank is as large as the model cardinality.
      is_model_identified = (rank_estimate >= previous_model_size);

      // Skip iteration if the model is not identified.
      if (!is_model_identified) continue;

    }

    // Compute the GCV score for the selected model
    score(par) = ComputeGCVScoreCpp(
      oracle_tfrp.col(par),
      idx_selected,
      covariance_factors_returns,
      variance_returns,
      mean_returns,
      scaling
    );

  }

  // Return the vector of GCV scores
  return score;

}

////////////////////////////////
///// CVTuningOracleTFRPCpp /////

arma::vec CVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const unsigned int n_folds
) {

  // Determine the number of observations
  const unsigned int n_observations = returns.n_rows;

  // Check the validity of n_folds
  if (n_folds == 0 || n_folds > n_observations) {
    Rcpp::stop("n_folds must be greater than 0 and less than or equal to the number of observations.");
  }

  // Calculate the fold size for cross-validation
  const unsigned int fold_size = std::max(1u, n_observations / n_folds);

  // Generate a vector of observation indices
  const arma::uvec observation_indices = arma::regspace<arma::uvec>(0, n_observations - 1);

  // Initialize a matrix to store prediction errors for each fold and penalty parameter
  arma::mat prediction_error(n_folds, penalty_parameters.n_elem);

  // Iterate over each fold for cross-validation
  for (unsigned int fold = 0; fold < n_folds; ++fold) {

    // Define the start and end indices for the test set in the current fold
    const unsigned int test_start = fold * fold_size;
    const unsigned int test_end = std::min(
      test_start + fold_size, n_observations
    ) - 1;

    // Determine the training indices (excluding the current test set)
    const arma::uvec train_indices = arma::find((observation_indices < test_start) || (observation_indices > test_end));

    // Extract training and test data for returns and factors
    const arma::mat& returns_train = returns.rows(train_indices);
    const arma::mat& factors_train = factors.rows(train_indices);
    const arma::mat& returns_test = returns.rows(test_start, test_end);

    // Compute covariance and mean of factors and returns for the training data
    const arma::mat cov_fac_ret_train = arma::cov(factors_train, returns_train);
    const arma::mat var_ret_train = arma::cov(returns_train);
    const arma::vec mean_ret_train = arma::mean(returns_train).t();

    // Compute adaptive weights for the training data
    const arma::vec weights_train = AdaptiveWeightsCpp(returns_train, factors_train, weighting_type);

    // Compute Oracle TFRP for the training data
    const arma::mat oracle_tfrp_train = OracleTFRPCpp(
      TFRPCpp(cov_fac_ret_train, var_ret_train, mean_ret_train),
      weights_train,
      penalty_parameters
    );

    // Iterate over each penalty parameter to compute prediction error
    for (unsigned int par = 0; par < penalty_parameters.n_elem; ++par) {

      prediction_error(fold, par) = ComputePredictionErrorCpp(
        oracle_tfrp_train.col(par),
        cov_fac_ret_train,
        var_ret_train,
        arma::mean(returns_test).t(),
        arma::var(returns_test).t()
      );

    }

  }

  // Return the mean prediction error across folds for each penalty parameter
  return arma::mean(prediction_error).t();

}

////////////////////////////////
///// RVTuningOracleTFRPCpp /////

arma::vec RVTuningOracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const unsigned int n_train_observations,
  const unsigned int n_test_observations,
  const unsigned int roll_shift
) {

  // Determine the total number of observations
  const unsigned int n_observations = returns.n_rows;

  // Check the validity of n_train_observations and n_test_observations
  if (n_train_observations == 0 || n_train_observations >= n_observations) {
    Rcpp::stop("n_train_observations must be greater than 0 and less than the total number of observations.");
  }
  if (n_test_observations == 0 || n_test_observations > n_observations / 2) {
    Rcpp::stop("n_test_observations must be greater than 0 and less than half the total number of observations.");
  }
  if (n_train_observations + n_test_observations > n_observations) {
    Rcpp::stop("Sum of n_train_observations and n_test_observations must not exceed the total number of observations.");
  }
  if (roll_shift == 0 || roll_shift > n_test_observations) {
    Rcpp::stop("roll_shift must be greater than 0 and less than the number of test observations.");
  }

  // Calculate the number of rolling windows
  const unsigned int n_rolls = (n_observations - n_train_observations - n_test_observations) / roll_shift + 1;

  // Initialize a matrix to store prediction errors for each roll and penalty parameter
  arma::mat prediction_error(n_rolls, penalty_parameters.n_elem);

  // Iterate over each roll for rolling validation
  for (unsigned int roll = 0; roll < n_rolls; ++roll) {

    // Define the start and end indices for the training and test sets in the current roll
    const unsigned int train_start = roll_shift * roll;
    const unsigned int train_end = n_train_observations - 1 + roll_shift * roll;
    const unsigned int test_start = n_train_observations + roll_shift * roll;
    const unsigned int test_end = roll < (n_rolls - 1) ?
    n_train_observations + n_test_observations - 1 + roll_shift * roll :
      n_observations - 1;

    // Extract training and test data for returns and factors
    const arma::mat returns_train = returns.rows(train_start, train_end);
    const arma::mat factors_train = factors.rows(train_start, train_end);
    const arma::mat returns_test = returns.rows(test_start, test_end);

    // Compute covariance and mean of factors and returns for the training data
    const arma::mat cov_fac_ret_train = arma::cov(factors_train, returns_train);
    const arma::mat var_ret_train = arma::cov(returns_train);
    const arma::vec mean_ret_train = arma::mean(returns_train).t();

    // Compute adaptive weights for the training data
    const arma::vec weights_train = AdaptiveWeightsCpp(returns_train, factors_train, weighting_type);

    // Compute Oracle TFRP for the training data
    const arma::mat oracle_tfrp_train = OracleTFRPCpp(
      TFRPCpp(cov_fac_ret_train, var_ret_train, mean_ret_train),
      weights_train,
      penalty_parameters
    );

    // Iterate over each penalty parameter to compute prediction error
    for (unsigned int par = 0; par < penalty_parameters.n_elem; ++par) {

      prediction_error(roll, par) = ComputePredictionErrorCpp(
        oracle_tfrp_train.col(par),
        cov_fac_ret_train,
        var_ret_train,
        arma::mean(returns_test).t(),
        arma::var(returns_test).t()
      );

    }

  }

  // Return the mean prediction error across rolls for each penalty parameter
  return arma::mean(prediction_error).t();

}

/////////////////////////////////////
///// ComputePredictionErrorCpp /////

double ComputePredictionErrorCpp(
  const arma::vec& oracle_tfrp_train,
  const arma::mat& covariance_factor_returns_train,
  const arma::mat& variance_returns_train,
  const arma::vec& mean_returns_test,
  const arma::vec& variance_returns_test
) {

  // Find the indices of the selected Oracle TFRP
  const arma::uvec idx_selected = arma::find(oracle_tfrp_train);

  // If no factors are selected, return the error assuming no model (only mean returns)
  if (idx_selected.empty()) {

    return arma::accu(arma::square(mean_returns_test) / variance_returns_test);

  }

  // Compute covariance of selected factors with returns
  const arma::mat cov_selected_fac_ret = covariance_factor_returns_train.rows(idx_selected);

  // Compute the inverse covariance of returns and covariance of returns with selected factors
  const arma::mat var_ret_inv_cov_ret_selected_fac = SolveSympd(
    variance_returns_train,
    cov_selected_fac_ret.t()
  );

  // Solve for beta coefficients for the selected factors
  const arma::mat beta_selected = SolveSympd(
    cov_selected_fac_ret * var_ret_inv_cov_ret_selected_fac,
    cov_selected_fac_ret
  ).t();

  // Calculate pricing error for the test data
  const arma::vec pricing_error_test = mean_returns_test - beta_selected * oracle_tfrp_train(idx_selected);

  // Return the weighted sum of squared pricing errors
  return arma::accu(arma::square(pricing_error_test) / variance_returns_test);

}

//////////////////////////////
///// ComputeGCVScoreCpp /////

double ComputeGCVScoreCpp(
  const arma::vec& oracle_tfrp,
  const arma::uvec& idx_selected,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const double scaling
) {

  // Compute covariance of selected factors with returns
  const arma::mat cov_selected_fac_ret = covariance_factors_returns.rows(idx_selected);

  const arma::mat var_ret_inv_cov_ret_selected_fac = SolveSympd(
    variance_returns,
    cov_selected_fac_ret.t()
  );

  // Solve for beta coefficients for the selected factors
  const arma::mat beta_selected = SolveSympd(
    cov_selected_fac_ret * var_ret_inv_cov_ret_selected_fac,
    cov_selected_fac_ret
  ).t();

  // Calculate pricing errors
  const arma::vec pricing_errors = mean_returns - beta_selected * oracle_tfrp(idx_selected);

  // Calculate the denominator of the GCV score formula
  const double denominator = std::max(
    .0001, 1. - static_cast<double>(idx_selected.n_elem) * scaling
  );

  // Return the GCV score which is the sum of squared pricing errors normalized by the square of the denominator
  return arma::accu(arma::square(pricing_errors)) / (denominator * denominator);

}
