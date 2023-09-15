// Author: Alberto Quaini

#include "tuning.h"
#include "tfrp.h"
#include "oracle_tfrp.h"
#include "adaptive_weights.h"
#include "identification_tests.h"

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

  // initialize the score to the inner product of the mean of returns
  // i.e., to the no model case
  arma::vec score(oracle_tfrp.n_cols, arma::fill::value(
    arma::dot(mean_returns, mean_returns)
  ));

  const double scaling = gcv_scaling_n_assets ?
  std::sqrt(mean_returns.n_elem) / n_observations :
    1. / n_observations;

  // in case the identification flag is true
  if (gcv_identification_check) {

    unsigned int cardinality_previous_model = factors.n_cols + 1;

    bool misidentification_flag = true;

    for (unsigned int par = 0; par < oracle_tfrp.n_cols; ++par) {

      const arma::uvec idx_selected = arma::find(oracle_tfrp.col(par));

       // if there are no selected factor, move to the next iteration
       if (!idx_selected.n_elem) continue;

       // check that the p-value of a beta rank test (Cheng Fang 2019) is above a
       // given threshold, indicating that we do not reject the null that the model
       // is not identified.
       // perform this check when: (a) the previously checked model was not
       // identified and (b) the model cardinality has increased
       if (misidentification_flag) {

         // if the model cardinality has not decreased, set the model score to
         // `score_no_model`
         if (idx_selected.n_elem >= cardinality_previous_model) continue;

         cardinality_previous_model = idx_selected.n_elem;

         const unsigned int kp2006_rank_estimate =
           IterativeKleibergenPaap2006BetaRankTestCpp(
             returns,
             factors.cols(idx_selected),
             target_level_kp2006_rank_test
           )["rank"];

         // if we do not reject the null that the model is not identified
         // set the model score to `score_no_model`
         // to decide, use a high threshold for the p-value
         if (kp2006_rank_estimate < factors.n_cols) continue;

         // otherwise set the mis-identification flag to false
         misidentification_flag = false;

        }

       score(par) = ComputeGCVScoreCpp(
         oracle_tfrp.col(par),
         idx_selected,
         covariance_factors_returns,
         variance_returns,
         mean_returns,
         scaling
       );

    }

  } else {

    for (unsigned int par = 0; par < oracle_tfrp.n_cols; ++par) {

      const arma::uvec idx_selected = arma::find(oracle_tfrp.col(par));

      // if there are no selected factor, move to the next iteration
      if (!idx_selected.n_elem) continue;

      score(par) = ComputeGCVScoreCpp(
        oracle_tfrp.col(par),
        idx_selected,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        scaling
      );

    }

  }

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

  const unsigned int n_observations = returns.n_rows;

  const unsigned int fold_size = n_observations / n_folds;

  const arma::uvec observation_indices =
    arma::regspace<arma::uvec>(0, n_observations - 1);

  arma::mat prediction_error(n_folds, penalty_parameters.n_elem);

  for (unsigned int fold = 0; fold < n_folds; ++fold) {

    const unsigned int test_start = fold * fold_size;

    const unsigned int test_end = fold < (n_folds - 1) ?
    test_start + fold_size - 1 :  n_observations - 1;

    const arma::uvec train_indices = arma::find(
      (observation_indices < test_start) or (observation_indices > test_end)
    );

    const arma::mat& returns_train = returns.rows(train_indices);
    const arma::mat& factors_train = factors.rows(train_indices);
    const arma::mat& returns_test = returns.rows(test_start, test_end);

    const arma::mat cov_fac_ret_train = arma::cov(factors_train, returns_train);
    const arma::mat var_ret_train = arma::cov(returns_train);
    const arma::vec mean_ret_train = arma::mean(returns_train).t();

    const arma::vec weights_train = AdaptiveWeightsCpp(
      returns_train,
      factors_train,
      weighting_type
    );

    const arma::mat oracle_tfrp_train = OracleTFRPCpp(
      TFRPCpp(cov_fac_ret_train, var_ret_train, mean_ret_train),
      weights_train,
      penalty_parameters
    );

    for (unsigned int par = 0; par < penalty_parameters.n_elem; ++par) {

      prediction_error(fold, par) += ComputePredictionErrorCpp(
        oracle_tfrp_train.col(par),
        cov_fac_ret_train,
        var_ret_train,
        arma::mean(returns_test).t(),
        arma::var(returns_test).t()
      );

    }

  }

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

  const unsigned int n_observations = returns.n_rows;

  const unsigned int n_rolls =
    (n_observations - n_train_observations - n_test_observations) / 12 + 1;

  arma::mat prediction_error(n_rolls, penalty_parameters.n_elem);

  for (unsigned int roll = 0; roll < n_rolls; ++roll) {

    const unsigned int train_start = roll_shift * roll;
    const unsigned int train_end =
      n_train_observations - 1 + roll_shift * roll;
    const unsigned int test_start =
      n_train_observations + roll_shift * roll;
    const unsigned int test_end = roll < (n_rolls - 1) ?
    n_train_observations + n_test_observations - 1 + roll_shift * roll :
      n_observations - 1;

    const arma::mat returns_train = returns.rows(train_start, train_end);
    const arma::mat factors_train = factors.rows(train_start, train_end);
    const arma::mat returns_test = returns.rows(test_start, test_end);

    const arma::mat cov_fac_ret_train = arma::cov(factors_train, returns_train);
    const arma::mat var_ret_train = arma::cov(returns_train);
    const arma::vec mean_ret_train = arma::mean(returns_train).t();

    const arma::vec weights_train = AdaptiveWeightsCpp(
      returns_train,
      factors_train,
      weighting_type
    );

    const arma::mat oracle_tfrp_train = OracleTFRPCpp(
      TFRPCpp(cov_fac_ret_train, var_ret_train, mean_ret_train),
      weights_train,
      penalty_parameters
    );

    for (unsigned int par = 0; par < penalty_parameters.n_elem; ++par) {

      prediction_error(roll, par) += ComputePredictionErrorCpp(
        oracle_tfrp_train.col(par),
        cov_fac_ret_train,
        var_ret_train,
        arma::mean(returns_test).t(),
        arma::var(returns_test).t()
      );

    }

  }

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

  const arma::uvec idx_selected = arma::find(oracle_tfrp_train);

  if (! idx_selected.n_elem) {

    return arma::dot(
      mean_returns_test, mean_returns_test / variance_returns_test
    );

  }

  const arma::mat cov_selected_fac_ret =
    covariance_factor_returns_train.rows(idx_selected);

  const arma::mat var_ret_inv_cov_ret_selected_fac = arma::solve(
    variance_returns_train, cov_selected_fac_ret.t(),
    arma::solve_opts::likely_sympd
  );

  const arma::mat beta_selected = arma::solve(
    cov_selected_fac_ret * var_ret_inv_cov_ret_selected_fac,
    cov_selected_fac_ret,
    arma::solve_opts::likely_sympd
  ).t();

  const arma::vec pricing_error_test = mean_returns_test -
    beta_selected * oracle_tfrp_train(idx_selected);

  return arma::dot(
    pricing_error_test,
    pricing_error_test / variance_returns_test
  );

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

  const arma::mat cov_selected_fac_ret =
    covariance_factors_returns.rows(idx_selected);

  const arma::mat var_ret_inv_cov_ret_selected_fac = arma::solve(
    variance_returns,
    cov_selected_fac_ret.t(),
    arma::solve_opts::likely_sympd
  );

  const arma::mat beta_selected = arma::solve(
    cov_selected_fac_ret * var_ret_inv_cov_ret_selected_fac,
    cov_selected_fac_ret,
    arma::solve_opts::likely_sympd
  ).t();

  const arma::vec pricing_errors = mean_returns - beta_selected *
    oracle_tfrp(idx_selected);

  const double denominator = 1. - std::min(
    (double)idx_selected.n_elem * scaling,
    .999
  );

  return arma::dot(pricing_errors, pricing_errors)
    / (denominator * denominator);

}
