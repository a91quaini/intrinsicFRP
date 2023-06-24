// Author: Alberto Quaini

#include "tuning_adaptive_ifrp.h"
#include "ifrp.h"
#include "adaptive_ifrp.h"
#include "adaptive_weights.h"
#include "identification_tests.h"

///////////////////////////////////
///// GCVScoreAdaptiveIFRPCpp /////

arma::vec GCVScoreAdaptiveIFRPCpp(
  const arma::mat& aifrp,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const unsigned int n_observations,
  const bool gcv_aic_scaling
) {

  const double score_no_model = arma::dot(mean_returns, mean_returns);

  arma::vec model_score(aifrp.n_cols, arma::fill::value(score_no_model));

  const double scaling = gcv_aic_scaling ? 1. / n_observations :
    std::log(n_observations) / n_observations;

  for (unsigned int par = 0; par < aifrp.n_cols; ++par) {

    const arma::uvec idx_selected = arma::find(aifrp.col(par));

    const unsigned int n_selected = idx_selected.n_elem;

    // if there are no selected factor, move to the next iteration
    if (! n_selected) continue;

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
      aifrp(idx_selected, arma::uvec(1, arma::fill::value(par)));

    const double denominator = 1. - std::min(idx_selected.n_elem * scaling, 1.);

    model_score(par) = arma::dot(pricing_errors, pricing_errors)
      / (denominator * denominator);

  }

  return model_score;

}

///////////////////////////////////
///// IGCVScoreAdaptiveIFRPCpp /////

arma::vec IGCVScoreAdaptiveIFRPCpp(
  const arma::mat& aifrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool gcv_aic_scaling
) {

  const double score_no_model = arma::dot(mean_returns, mean_returns);

  arma::vec model_score(aifrp.n_cols, arma::fill::value(score_no_model));

  const unsigned int n_observations = factors.n_rows;

  const double scaling = gcv_aic_scaling ? 1. / n_observations :
    std::log(n_observations) / n_observations;

  unsigned int cardinality_previous_model = factors.n_cols + 1;

  bool misidentification_flag = true;

  for (unsigned int par = 0; par < aifrp.n_cols; ++par) {

    const arma::uvec idx_selected = arma::find(aifrp.col(par));

    // if there are no selected factor, move to the next iteration
    if (! idx_selected.n_elem) continue;

    const arma::mat cov_selected_fac_ret =
      covariance_factors_returns.rows(idx_selected);

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

      const arma::mat beta_selected = arma::solve(
        arma::cov(factors.cols(idx_selected)),
        cov_selected_fac_ret,
        arma::solve_opts::likely_sympd
      ).t();

      const arma::vec2 identification_statistics =
        BetaRankChenFang2019StatisticAndPvalueCpp(
          returns,
          factors.cols(idx_selected),
          beta_selected
        );

      Rcpp::Rcout << "p-value = " << identification_statistics(1) <<"\n";

      // if we do not reject the null that the model is not identified
      // set the model score to `score_no_model`
      // to decide, use a high threshold for the p-value
      if (identification_statistics(1) >= 0.1) continue;

      // otherwise set the mis-identification flag to false
      misidentification_flag = false;

    }

    const arma::mat var_ret_inv_cov_ret_selected_fac = arma::solve(
      variance_returns,
      cov_selected_fac_ret.t(),
      arma::solve_opts::likely_sympd
    );

    const arma::mat beta_intrinsic = arma::solve(
      cov_selected_fac_ret * var_ret_inv_cov_ret_selected_fac,
      cov_selected_fac_ret,
      arma::solve_opts::likely_sympd
    ).t();

    const arma::vec pricing_errors = mean_returns - beta_intrinsic *
      aifrp(idx_selected, arma::uvec(1, arma::fill::value(par)));

    const double denominator = 1. - std::min(idx_selected.n_elem * scaling, 1.);

    model_score(par) = arma::dot(pricing_errors, pricing_errors)
      / (denominator * denominator);

  }

  return model_score;

}

///////////////////////////////////////////
///// WeightedGCVScoreAdaptiveIFRPCpp /////

arma::vec WeightedGCVScoreAdaptiveIFRPCpp(
  const arma::mat& aifrp,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const unsigned int n_observations,
  const bool gcv_aic_scaling
) {

  const double score_no_model = arma::dot(mean_returns, arma::solve(
    variance_returns, mean_returns, arma::solve_opts::likely_sympd
  ));

  arma::vec model_score(aifrp.n_cols, arma::fill::value(score_no_model));

  const double scaling = gcv_aic_scaling ? 1. / n_observations :
    std::log(n_observations) / n_observations;

  for (unsigned int par = 0; par < aifrp.n_cols; ++par) {

    const arma::uvec idx_selected = arma::find(aifrp.col(par));

    const unsigned int n_selected = idx_selected.n_elem;

    // if there are no selected factor, move to the next iteration
    if (! n_selected) continue;

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
      aifrp(idx_selected, arma::uvec(1, arma::fill::value(par)));

    const double denominator = 1. - std::min(idx_selected.n_elem * scaling, 1.);

    model_score(par) = arma::dot(pricing_errors, arma::solve(
      variance_returns, pricing_errors, arma::solve_opts::likely_sympd
    )) / (denominator * denominator);

  }

  return model_score;

}

//////////////////////////////////
///// CVScoreAdaptiveIFRPCpp /////

arma::vec CVScoreAdaptiveIFRPCpp(
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

  arma::vec model_score(penalty_parameters.n_elem);

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

    const arma::mat aifrp_train = AdaptiveIFRPCpp(
      IFRPCpp(cov_fac_ret_train, var_ret_train, mean_ret_train),
      weights_train,
      penalty_parameters
    );

    for (unsigned int par = 0; par < penalty_parameters.n_elem; ++par) {

      model_score(par) += ComputePredictionErrorCpp(
        aifrp_train.col(par),
        cov_fac_ret_train,
        var_ret_train,
        arma::mean(returns_test).t(),
        arma::var(returns_test).t()
      );

    }

  }

  return model_score;

}

//////////////////////////////////
///// RVScoreAdaptiveIFRPCpp /////

arma::vec RVScoreAdaptiveIFRPCpp(
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

  arma::vec model_score(penalty_parameters.n_elem);

  for (unsigned int roll=0; roll < n_rolls; ++roll) {

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

    const arma::mat aifrp_train = AdaptiveIFRPCpp(
      IFRPCpp(cov_fac_ret_train, var_ret_train, mean_ret_train),
      weights_train,
      penalty_parameters
    );

    for (unsigned int par = 0; par < penalty_parameters.n_elem; ++par) {

      model_score(par) += ComputePredictionErrorCpp(
        aifrp_train.col(par),
        cov_fac_ret_train,
        var_ret_train,
        arma::mean(returns_test).t(),
        arma::var(returns_test).t()
      );

    }

  }

  return model_score;

}

/////////////////////////////////////
///// ComputePredictionErrorCpp /////

double ComputePredictionErrorCpp(
  const arma::vec& aifrp_train,
  const arma::mat& covariance_factor_returns_train,
  const arma::mat& variance_returns_train,
  const arma::vec& mean_returns_test,
  const arma::vec& variance_returns_test
) {

  const arma::uvec idx_selected = arma::find(aifrp_train);

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
    beta_selected * aifrp_train(idx_selected);

  return arma::dot(
    pricing_error_test,
    pricing_error_test / variance_returns_test
  );

}
