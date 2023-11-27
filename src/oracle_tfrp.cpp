// Author: Alberto Quaini

#include "oracle_tfrp.h"
#include "tfrp.h"
#include "adaptive_weights.h"
#include "tuning.h"

Rcpp::List OracleTFRPCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::vec& penalty_parameters,
  const char weighting_type,
  const char tuning_type,
  const bool one_stddev_rule,
  const bool gcv_scaling_n_assets,
  const bool gcv_identification_check,
  const double target_level_kp2006_rank_test,
  const unsigned int n_folds,
  const unsigned int n_train_observations,
  const unsigned int n_test_observations,
  const unsigned int roll_shift,
  const bool relaxed,
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  const Rcpp::List output = [&] {
    // Switch case for different tuning types
    switch (tuning_type) {
    case 'g':
     return OracleTFRPGCVCpp(
        returns,
        factors,
        arma::cov(factors, returns),
        arma::cov(returns),
        arma::mean(returns).t(),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        gcv_scaling_n_assets,
        gcv_identification_check,
        target_level_kp2006_rank_test,
        relaxed,
        include_standard_errors,
        hac_prewhite
      );
    case 'c':
      return OracleTFRPCVCpp(
        returns,
        factors,
        arma::cov(factors, returns),
        arma::cov(returns),
        arma::mean(returns).t(),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        n_folds,
        relaxed,
        include_standard_errors,
        hac_prewhite
      );
    case 'r':
      return OracleTFRPRVCpp(
        returns,
        factors,
        arma::cov(factors, returns),
        arma::cov(returns),
        arma::mean(returns).t(),
        penalty_parameters,
        weighting_type,
        one_stddev_rule,
        n_train_observations,
        n_test_observations,
        roll_shift,
        relaxed,
        include_standard_errors,
        hac_prewhite
      );
    default:
      Rcpp::stop("Invalid tuning type");
    }
  }();

  // Return output list.
  return output;

}

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
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  // Compute the Oracle TFRP for a range of penalty parameters
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

  // Compute the GCV scores for the Oracle TFRP
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

  // Identify the optimal penalty parameter based on GCV score
  const unsigned int idx_optimal_parameter = one_stddev_rule ?
    ComputeOneStdDevRuleCpp(score) : score.index_min();

  // Return the Oracle TFRP and standard errors if requested
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
    mean_returns,
    hac_prewhite
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
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  // Compute CV scores for Oracle TFRP
  const arma::vec score = CVTuningOracleTFRPCpp(
    returns,
    factors,
    penalty_parameters,
    weighting_type,
    n_folds
  );

  // Identify the optimal penalty parameter based on CV score
  const unsigned int idx_optimal_parameter = one_stddev_rule ?
  ComputeOneStdDevRuleCpp(score) : score.index_min();

  // Compute the Oracle TFRP for the optimal penalty parameter
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

  // Return the Oracle TFRP and standard errors if requested
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
    mean_returns,
    hac_prewhite
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
  const bool include_standard_errors,
  const bool hac_prewhite
) {

  // Compute RV scores for Oracle TFRP
  const arma::vec score = RVTuningOracleTFRPCpp(
    returns,
    factors,
    penalty_parameters,
    weighting_type,
    n_train_observations,
    n_test_observations,
    roll_shift
  );

  // Identify the optimal penalty parameter based on RV score
  const unsigned int idx_optimal_parameter = one_stddev_rule ?
  ComputeOneStdDevRuleCpp(score) : score.index_min();

  // Compute the Oracle TFRP for the optimal penalty parameter
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

  // Return the Oracle TFRP and standard errors if requested
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
    mean_returns,
    hac_prewhite
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

  // Calculate sign of the tradable factor risk premia
  const arma::vec sign_factor_rp = arma::sign(tradable_frp);
  // Initialize a temporary matrix for computing Oracle TFRP
  arma::mat temp = - weights * penalty_parameters.t();
  // Add the absolute value of tradable factor risk premia to each column of temp
  temp.each_col() += sign_factor_rp % tradable_frp;
  // Clamp the values in temp to non-negative values
  temp.clamp(0., arma::datum::inf);

  // Return the sign-adjusted Oracle TFRP
  return sign_factor_rp % temp.each_col();

}

arma::vec OracleTFRPCpp(
  const arma::vec& tradable_frp,
  const arma::vec& weights,
  const double penalty_parameter
) {

  // Calculate sign of the tradable factor risk premia
  const arma::vec sign_ifrp = arma::sign(tradable_frp);

  // Compute the Oracle TFRP by soft-thresholding the tradable factor risk premia
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

  // Initialize a vector for relaxed Oracle TFRP with the number of rows equal to covariance_factors_returns
  arma::vec relaxed_otfrp(covariance_factors_returns.n_rows);

  // Compute TFRP only for the indices that are non-zero in the Oracle TFRP
  relaxed_otfrp(idx_nonzero) = TFRPCpp(
    covariance_factors_returns.rows(idx_nonzero),
    variance_returns,
    mean_returns
  );

  // Return the relaxed Oracle TFRP
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
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Return zeros if no non-zero indices are present (no premia to compute errors for)
  if (!idx_nonzero.n_elem) { return arma::zeros(factors.n_cols); }

  // Initialize a vector for standard errors
  arma::vec standard_errors(factors.n_cols);

  // Compute standard errors only for the non-zero indices in the Oracle TFRP
  standard_errors(idx_nonzero) = StandardErrorsTFRPCpp(
    returns,
    factors.cols(idx_nonzero),
    covariance_factors_returns.rows(idx_nonzero),
    variance_returns,
    mean_returns,
    hac_prewhite
  );

  // Return the standard errors for Oracle TFRP
  return standard_errors;

}

///////////////////////////////////
///// ComputeOneStdDevRuleCpp /////

unsigned int ComputeOneStdDevRuleCpp(
  const arma::vec& score
) {

  // Find the index of the minimum score
  unsigned int idx_optimal_parameter = score.index_min();

  // Consider the scores to the right of the minimum score
  const arma::vec score_right_of_min = score(
    arma::span(idx_optimal_parameter, score.n_elem - 1)
  );

  // Apply the one standard deviation rule to find the most parsimonious model
  // whose score is not more than one standard deviation above the minimum score
  return idx_optimal_parameter + arma::max(arma::find(
      score_right_of_min <= arma::min(score_right_of_min) +
        arma::stddev(score_right_of_min)
  ));

}
