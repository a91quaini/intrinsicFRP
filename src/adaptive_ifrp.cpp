// Author: Alberto Quaini

#include "adaptive_ifrp.h"
#include "ifrp.h"
#include "adaptive_weights.h"
#include "tuning_adaptive_ifrp.h"
#include "hac_standard_errors.h"

////////////////////////////////////
///// OptimalAdaptiveIFRPGCVCpp ////

Rcpp::List OptimalAdaptiveIFRPGCVCpp(
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
  const double level_kp2006_rank_test
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

  return Rcpp::List::create(
    Rcpp::Named("risk_premia") = aifrp.col(idx_optimal_parameter),
    Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
    Rcpp::Named("model_score") = model_score
  );

}

////////////////////////////////
///// OptimalAdaptiveIFRPCV ////

Rcpp::List OptimalAdaptiveIFRPCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const bool one_stddev_rule,
  const unsigned int n_folds
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

  return Rcpp::List::create(
    Rcpp::Named("risk_premia") = AdaptiveIFRPCpp(
      IFRPCpp(
        arma::cov(factors, returns),
        arma::cov(returns),
        arma::mean(returns).t()
      ),
      AdaptiveWeightsCpp(
        returns,
        factors,
        weighting_type
      ),
      penalty_parameters(idx_optimal_parameter)
    ),
    Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
    Rcpp::Named("model_score") = model_score
  );

}

///////////////////////////////////
///// OptimalAdaptiveIFRPRVCpp ////

Rcpp::List OptimalAdaptiveIFRPRVCpp(
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
  const unsigned int roll_shift
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

  return Rcpp::List::create(
    Rcpp::Named("risk_premia") = AdaptiveIFRPCpp(
      IFRPCpp(
        arma::cov(factors, returns),
        arma::cov(returns),
        arma::mean(returns).t()
      ),
      AdaptiveWeightsCpp(
        returns,
        factors,
        weighting_type
      ),
      penalty_parameters(idx_optimal_parameter)
    ),
    Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
    Rcpp::Named("model_score") = model_score
  );

}

///////////////////////////
///// AdaptiveIFRPCpp /////

arma::mat AdaptiveIFRPCpp(
  const arma::vec& ifrp,
  const arma::vec& weights,
  const arma::vec& penalty_parameters
) {

  const arma::vec sign_factor_rp = arma::sign(ifrp);
  arma::mat temp = - weights * penalty_parameters.t();
  temp.each_col() += sign_factor_rp % ifrp;
  temp.clamp(0., arma::datum::inf);

  return sign_factor_rp % temp.each_col();

}

arma::vec AdaptiveIFRPCpp(
  const arma::vec& ifrp,
  const arma::vec& weights,
  const double penalty_parameter
) {

  const arma::vec sign_ifrp = arma::sign(ifrp);

  return sign_ifrp % arma::clamp(
    sign_ifrp % ifrp - penalty_parameter * weights,
    0., arma::datum::inf
  );

}

/////////////////////////////////////////
///// StandardErrorsAdaptiveIFRPCpp /////

arma::vec StandardErrorsAdaptiveIFRPCpp(
  const arma::vec& aifrp,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& mean_factors
) {

  const arma::uvec idx_selected_rp = arma::find(aifrp);

  if (! idx_selected_rp.n_elem) {return arma::zeros(mean_factors.n_elem);}

  const arma::mat var_ret_inv_cov_ret_selected_fac = arma::solve(
    variance_returns,
    covariance_factors_returns.rows(idx_selected_rp).t(),
    arma::solve_opts::likely_sympd
  );

  const arma::mat returns_centred = returns.each_row() - mean_returns.t();

  const arma::mat selected_factors = factors.cols(idx_selected_rp);

  const arma::mat selected_factors_centred =
    selected_factors.each_row() - mean_factors(idx_selected_rp).t();

  const arma::vec ret_cen_var_ret_inv_mean_ret = returns_centred * arma::solve(
    variance_returns,
    mean_returns,
    arma::solve_opts::likely_sympd
  );

  const arma::mat ret_cen_var_ret_inv_cov_ret_selected_fac =
    returns_centred * var_ret_inv_cov_ret_selected_fac;

  arma::vec standard_errors(factors.n_cols);

  standard_errors(idx_selected_rp) = HACStandardErrorsCpp(
    // covariance term
    selected_factors_centred.each_col() % ret_cen_var_ret_inv_mean_ret -
    // variance term
    ret_cen_var_ret_inv_cov_ret_selected_fac.each_col() %
    ret_cen_var_ret_inv_mean_ret +
    // mean term
    returns_centred * var_ret_inv_cov_ret_selected_fac
  );

  return standard_errors;

}
