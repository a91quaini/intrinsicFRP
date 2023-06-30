// Author: Alberto Quaini

#include "adaptive_frp.h"
#include "adaptive_ifrp.h"
#include "ifrp.h"
#include "frp.h"
#include "adaptive_weights.h"
#include "tuning_adaptive_ifrp.h"
#include "hac_standard_errors.h"

////////////////////////////////////
///// OptimalAdaptiveFRPGCVCpp ////

Rcpp::List OptimalAdaptiveFRPGCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const bool one_stddev_rule,
  const bool gcv_vr_weighting,
  const bool gcv_scaling_n_assets,
  const bool gcv_identification_check,
  const unsigned int n_bootstrap_cf2019_rank_test,
  const double level_cf2019_rank_test,
  const double level_kp2006_rank_test,
  const bool include_standard_errors
) {

  arma::mat aifrp = AdaptiveIFRPCpp(
    IFRPCpp(
      covariance_factors_returns,
      variance_returns,
      mean_returns
    ),
    AdaptiveWeightsCpp(
      returns,
      factors,
      weighting_type
    ),
    penalty_parameters
  );

  arma::vec model_score(penalty_parameters.n_elem);

  if (gcv_vr_weighting) {

    model_score = gcv_identification_check ?
      WeightedIGCVScoreAdaptiveIFRPCpp(
        aifrp,
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        returns.n_rows,
        gcv_scaling_n_assets,
        n_bootstrap_cf2019_rank_test,
        level_cf2019_rank_test,
        level_kp2006_rank_test
      ) :
      WeightedGCVScoreAdaptiveIFRPCpp(
        aifrp,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        returns.n_rows,
        gcv_scaling_n_assets
      );

  } else {

    model_score = gcv_identification_check ?
      IGCVScoreAdaptiveIFRPCpp(
        aifrp,
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        gcv_scaling_n_assets,
        n_bootstrap_cf2019_rank_test,
        level_cf2019_rank_test,
        level_kp2006_rank_test
      ) :
      GCVScoreAdaptiveIFRPCpp(
        aifrp,
        covariance_factors_returns,
        variance_returns,
        mean_returns,
        returns.n_rows,
        gcv_scaling_n_assets
      );

  }

  unsigned int idx_optimal_parameter = model_score.index_min();

  if (one_stddev_rule) {

    const arma::vec model_score_right_of_min = model_score(
      arma::span(idx_optimal_parameter, model_score.n_elem - 1)
    );

    idx_optimal_parameter += arma::max(arma::find(
      model_score_right_of_min <= arma::min(model_score_right_of_min) +
      arma::stddev(model_score_right_of_min)
    ));

  }

  const arma::uvec idx_selected = arma::find(aifrp.col(idx_optimal_parameter));

  arma::vec afrp(factors.n_cols);
  afrp(idx_selected) = KRSFRPCpp(
    arma::solve(
      arma::cov(factors.cols(idx_selected)),
      covariance_factors_returns.rows(idx_selected),
      arma::solve_opts::likely_sympd
    ).t(),
    mean_returns,
    variance_returns
  );

  return include_standard_errors ?
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = afrp,
      Rcpp::Named("standard_errors") = StandardErrorsAdaptiveFRPCpp(
        afrp,
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = model_score
    ) :
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = afrp,
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = model_score
    );

}

////////////////////////////////
///// OptimalAdaptiveFRPCV ////

Rcpp::List OptimalAdaptiveFRPCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const bool one_stddev_rule,
  const unsigned int n_folds,
  const bool include_standard_errors
) {

  const arma::vec model_score = CVScoreAdaptiveIFRPCpp(
    returns,
    factors,
    penalty_parameters,
    weighting_type,
    n_folds
  );

  unsigned int idx_optimal_parameter = model_score.index_min();

  if (one_stddev_rule) {

    const arma::vec model_score_right_of_min = model_score(
      arma::span(idx_optimal_parameter, model_score.n_elem - 1)
    );

    idx_optimal_parameter += arma::max(arma::find(
      model_score_right_of_min <= arma::min(model_score_right_of_min) +
      arma::stddev(model_score_right_of_min)
    ));

  }

  const arma::vec aifrp = AdaptiveIFRPCpp(
    IFRPCpp(
      covariance_factors_returns,
      variance_returns,
      mean_returns
    ),
    AdaptiveWeightsCpp(
      returns,
      factors,
      weighting_type
    ),
    penalty_parameters(idx_optimal_parameter)
  );

  const arma::uvec idx_selected = arma::find(aifrp.col(idx_optimal_parameter));

  arma::vec afrp(factors.n_cols);
  afrp(idx_selected) = KRSFRPCpp(
    arma::solve(
      arma::cov(factors.cols(idx_selected)),
      covariance_factors_returns.rows(idx_selected),
      arma::solve_opts::likely_sympd
    ).t(),
    mean_returns,
    variance_returns
  );

  return include_standard_errors ?
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = afrp,
      Rcpp::Named("standard_errors") = StandardErrorsAdaptiveFRPCpp(
        afrp,
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = model_score
    ) :
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = afrp,
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = model_score
    );

}

///////////////////////////////////
///// OptimalAdaptiveFRPRVCpp ////

Rcpp::List OptimalAdaptiveFRPRVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const bool one_stddev_rule,
  const unsigned int n_train_observations,
  const unsigned int n_test_observations,
  const unsigned int roll_shift,
  const bool include_standard_errors
) {

  const arma::vec model_score = RVScoreAdaptiveIFRPCpp(
    returns,
    factors,
    penalty_parameters,
    weighting_type,
    n_train_observations,
    n_test_observations,
    roll_shift
  );

  unsigned int idx_optimal_parameter = model_score.index_min();

  if (one_stddev_rule) {

    const arma::vec model_score_right_of_min = model_score(
      arma::span(idx_optimal_parameter, model_score.n_elem - 1)
    );

    idx_optimal_parameter += arma::max(arma::find(
      model_score_right_of_min <= arma::min(model_score_right_of_min) +
      arma::stddev(model_score_right_of_min)
    ));

  }

  const arma::vec aifrp = AdaptiveIFRPCpp(
    IFRPCpp(
      covariance_factors_returns,
      variance_returns,
      mean_returns
    ),
    AdaptiveWeightsCpp(
      returns,
      factors,
      weighting_type
    ),
    penalty_parameters(idx_optimal_parameter)
  );

  const arma::uvec idx_selected = arma::find(aifrp.col(idx_optimal_parameter));

  arma::vec afrp(factors.n_cols);
  afrp(idx_selected) = KRSFRPCpp(
    arma::solve(
      arma::cov(factors.cols(idx_selected)),
      covariance_factors_returns.rows(idx_selected),
      arma::solve_opts::likely_sympd
    ).t(),
    mean_returns,
    variance_returns
  );

  return include_standard_errors ?
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = afrp,
      Rcpp::Named("standard_errors") = StandardErrorsAdaptiveFRPCpp(
        afrp,
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = model_score
    ) :
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = afrp,
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = model_score
    );

}

/////////////////////////////////////////
///// StandardErrorsAdaptiveFRPCpp /////

arma::vec StandardErrorsAdaptiveFRPCpp(
  const arma::vec& afrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  const arma::uvec idx_selected = arma::find(afrp);

  if (! idx_selected.n_elem) {return arma::zeros(factors.n_cols);}

  const arma::mat factors_selected = factors.cols(idx_selected);

  arma::vec standard_errors(factors.n_cols);
  standard_errors(idx_selected) = StandardErrorsKRSFRPCpp(
    afrp(idx_selected),
    returns,
    factors_selected,
    arma::solve(
      arma::cov(factors_selected),
      covariance_factors_returns.rows(idx_selected),
      arma::solve_opts::likely_sympd
    ).t(),
    covariance_factors_returns.rows(idx_selected),
    variance_returns,
    mean_returns
  );

  return standard_errors;

}
