// Author: Alberto Quaini

#include "oracle_tfrp.h"
#include "tfrp.h"
#include "adaptive_weights.h"
#include "tuning.h"
#include "hac_standard_errors.h"

////////////////////////////////////
///// OracleTFRPGCVCpp ////

Rcpp::List OracleTFRPGCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const bool one_stddev_rule,
  const bool gcv_scaling_n_assets,
  const bool gcv_identification_check,
  const double target_level_kp2006_rank_test,
  const bool relaxed,
  const bool include_standard_errors
) {

  const arma::mat oracle_tfrp = OracleTFRPCpp(
    TFRPCpp(
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

  const arma::vec score = GCVTuningOracleTFRPCpp(
    returns,
    factors,
    oracle_tfrp,
    covariance_factors_returns,
    variance_returns,
    mean_returns,
    returns.n_rows,
    gcv_scaling_n_assets,
    gcv_identification_check,
    target_level_kp2006_rank_test
  );

  const unsigned int idx_optimal_parameter = one_stddev_rule ?
    ComputeOneStdDevRuleCpp(score) : score.index_min();

  return include_standard_errors ?
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = relaxed ?
        RelaxOracleTFRP(
          arma::find(oracle_tfrp.col(idx_optimal_parameter)),
          covariance_factors_returns,
          variance_returns,
          mean_returns
        ) : oracle_tfrp.col(idx_optimal_parameter),
      Rcpp::Named("standard_errors") = StandardErrorsOracleTFRPCpp(
        arma::find(oracle_tfrp.col(idx_optimal_parameter)),
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = score
    ) :
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = relaxed ?
        RelaxOracleTFRP(
          arma::find(oracle_tfrp.col(idx_optimal_parameter)),
          covariance_factors_returns,
          variance_returns,
          mean_returns
        ) : oracle_tfrp.col(idx_optimal_parameter),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = score
    );

}

///////////////////////
///// OracleTFRPCV ////

Rcpp::List OracleTFRPCVCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const bool one_stddev_rule,
  const unsigned int n_folds,
  const bool relaxed,
  const bool include_standard_errors
) {

  const arma::vec score = CVTuningOracleTFRPCpp(
    returns,
    factors,
    penalty_parameters,
    weighting_type,
    n_folds
  );

  const unsigned int idx_optimal_parameter = one_stddev_rule ?
    ComputeOneStdDevRuleCpp(score) : score.index_min();

  arma::vec oracle_tfrp = OracleTFRPCpp(
    TFRPCpp(
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

  return include_standard_errors ?
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = relaxed ?
        RelaxOracleTFRP(
          arma::find(oracle_tfrp),
          covariance_factors_returns,
          variance_returns,
          mean_returns
        ) : oracle_tfrp,
      Rcpp::Named("standard_errors") = StandardErrorsOracleTFRPCpp(
        arma::find(oracle_tfrp),
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = score
    ) :
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = relaxed ?
        RelaxOracleTFRP(
          arma::find(oracle_tfrp),
          covariance_factors_returns,
          variance_returns,
          mean_returns
        ) : oracle_tfrp,
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = score
    );

}

//////////////////////////
///// OracleTFRPRVCpp ////

Rcpp::List OracleTFRPRVCpp(
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
  const bool relaxed,
  const bool include_standard_errors
) {

  const arma::vec score = RVTuningOracleTFRPCpp(
    returns,
    factors,
    penalty_parameters,
    weighting_type,
    n_train_observations,
    n_test_observations,
    roll_shift
  );

  const unsigned int idx_optimal_parameter = one_stddev_rule ?
    ComputeOneStdDevRuleCpp(score) : score.index_min();

  const arma::vec oracle_tfrp = OracleTFRPCpp(
    TFRPCpp(
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

  return include_standard_errors ?
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = relaxed ?
        RelaxOracleTFRP(
          arma::find(oracle_tfrp),
          covariance_factors_returns,
          variance_returns,
          mean_returns
        ) : oracle_tfrp,
      Rcpp::Named("standard_errors") = StandardErrorsOracleTFRPCpp(
        arma::find(oracle_tfrp),
        returns,
        factors,
        covariance_factors_returns,
        variance_returns,
        mean_returns
      ),
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = score
    ) :
    Rcpp::List::create(
      Rcpp::Named("risk_premia") = relaxed ?
        RelaxOracleTFRP(
          arma::find(oracle_tfrp),
          covariance_factors_returns,
          variance_returns,
          mean_returns
        ) : oracle_tfrp,
      Rcpp::Named("penalty_parameter") = penalty_parameters(idx_optimal_parameter),
      Rcpp::Named("model_score") = score
    );

}

/////////////////////////
///// OracleTFRPCpp /////

arma::mat OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const arma::vec& penalty_parameters
) {

  const arma::vec sign_factor_rp = arma::sign(tradable_frp);
  arma::mat temp = - weights * penalty_parameters.t();
  temp.each_col() += sign_factor_rp % tradable_frp;
  temp.clamp(0., arma::datum::inf);

  return sign_factor_rp % temp.each_col();

}

arma::vec OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const double penalty_parameter
) {

  const arma::vec sign_ifrp = arma::sign(tradable_frp);

  return sign_ifrp % arma::clamp(
    sign_ifrp % tradable_frp - penalty_parameter * weights,
    0., arma::datum::inf
  );

}

//////////////////////////////
///// RelaxOracleTFRPCpp /////

arma::vec RelaxOracleTFRP(
  const arma::uvec& idx_nonzero,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  arma::vec relaxed_otfrp(covariance_factors_returns.n_rows);

  relaxed_otfrp(idx_nonzero) = TFRPCpp(
    covariance_factors_returns.rows(idx_nonzero),
    variance_returns,
    mean_returns
  );

  return relaxed_otfrp;

}

/////////////////////////////////////////
///// StandardErrorsOracleTFRPCpp /////

arma::vec StandardErrorsOracleTFRPCpp(
  const arma::uvec& idx_nonzero,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_factors_returns,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  if (!idx_nonzero.n_elem) {return arma::zeros(factors.n_cols);}

  arma::vec standard_errors(factors.n_cols);

  standard_errors(idx_nonzero) = StandardErrorsTFRPCpp(
    returns,
    factors.cols(idx_nonzero),
    covariance_factors_returns.rows(idx_nonzero),
    variance_returns,
    mean_returns
  );

  return standard_errors;

}

///////////////////////////////////
///// ComputeOneStdDevRuleCpp /////

unsigned int ComputeOneStdDevRuleCpp(
  const arma::vec& score
) {

    unsigned int idx_optimal_parameter = score.index_min();

    const arma::vec score_right_of_min = score(
      arma::span(idx_optimal_parameter, score.n_elem - 1)
    );

    return idx_optimal_parameter + arma::max(arma::find(
      score_right_of_min <= arma::min(score_right_of_min) +
        arma::stddev(score_right_of_min)
    ));

}
