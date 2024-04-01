// Author: Alberto Quaini

#include "sdf_coefficients.h"
#include "hac_covariance.h"
#include "gkr_factor_screening.h"
#include "utils.h"

/////////////////////////////
///// SDFCoefficientsCpp ////

Rcpp::List SDFCoefficientsCpp(
  const arma::mat& returns,
  const arma::mat& factors,
  const bool misspecification_robust,
  const bool include_standard_errors,
  const bool hac_prewhite,
  const double target_level_gkr2014_screening
) {

  // If target_level_gkr2014_screening > 0, perform the GKR screening procedure
  if (target_level_gkr2014_screening > 0.) {

    // Store the indices of the selected factors
    const arma::uvec selected_factor_indices = static_cast<arma::uvec>(
      GKRFactorScreeningCpp(
        returns,
        factors,
        target_level_gkr2014_screening,
        hac_prewhite
      )["selected_factor_indices"]
    );

    // If the GKR procedure removed all the factors, return a list containing empty vectors.
    if (selected_factor_indices.empty()) return include_standard_errors ?
      Rcpp::List::create(
        Rcpp::Named("sdf_coefficients") = arma::vec(),
        Rcpp::Named("standard_errors") = arma::vec(),
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      ) :
      Rcpp::List::create(
        Rcpp::Named("sdf_coefficients") = arma::vec(),
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      );

    // Otherwise, compute the SDF coefficients and eventual standard errors
    const Rcpp::List output = ReturnSDFCoefficientsCpp(
      returns,
      factors.cols(selected_factor_indices),
      misspecification_robust,
      include_standard_errors,
      hac_prewhite
    );

    // Return the output list with additionally the selected factor indices
    return include_standard_errors ?
      Rcpp::List::create(
        Rcpp::Named("sdf_coefficients") = output["sdf_coefficients"],
        Rcpp::Named("standard_errors") = output["standard_errors"],
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      ) :
      Rcpp::List::create(
        Rcpp::Named("sdf_coefficients") = output["sdf_coefficients"],
        Rcpp::Named("selected_factor_indices") = selected_factor_indices
      );

  }

  // Otherwise, return list of results without selected factor indices
  return ReturnSDFCoefficientsCpp(
    returns,
    factors,
    misspecification_robust,
    include_standard_errors,
    hac_prewhite
  );

}

///////////////////////////////////
///// ReturnSDFCoefficientsCpp ////

Rcpp::List ReturnSDFCoefficientsCpp(
    const arma::mat& returns,
    const arma::mat& factors,
    const bool misspecification_robust,
    const bool include_standard_errors,
    const bool hac_prewhite
) {

  // Check if standard errors are to be included in the calculation
  if (include_standard_errors) {

    // Compute covariance matrices and mean returns
    const arma::mat covariance_returns_factors = arma::cov(returns, factors);
    const arma::mat variance_returns = arma::cov(returns);
    const arma::vec mean_returns = arma::mean(returns).t();

    // Compute SDF coefficients based on whether misspecification is robust or not
    const arma::vec sdf_coefficients = misspecification_robust ?
      GKRSDFCoefficientsCpp(covariance_returns_factors, variance_returns, mean_returns) :
      FMSDFCoefficientsCpp(covariance_returns_factors, mean_returns);

    // Return SDF coefficients and standard errors in a list
    return Rcpp::List::create(
      Rcpp::Named("sdf_coefficients") =  sdf_coefficients,
      Rcpp::Named("standard_errors") = misspecification_robust ?
        StandardErrorsGKRSDFCoefficientsCpp(
          sdf_coefficients,
          returns,
          factors,
          covariance_returns_factors,
          variance_returns,
          mean_returns,
          hac_prewhite
        ) :
        StandardErrorsFMSDFCoefficientsCpp(
          sdf_coefficients,
          returns,
          factors,
          covariance_returns_factors,
          mean_returns,
          hac_prewhite
        )
    );

  }

    // Return SDF coefficients and standard errors in a list
    return Rcpp::List::create(
      Rcpp::Named("sdf_coefficients") =  misspecification_robust ?
        GKRSDFCoefficientsCpp(arma::cov(returns, factors),  arma::cov(returns), arma::mean(returns).t()) :
        FMSDFCoefficientsCpp(arma::cov(returns, factors), arma::mean(returns).t())
    );

}

///////////////////////////////
///// FMSDFCoefficientsCpp ////

arma::vec FMSDFCoefficientsCpp(
    const arma::mat& covariance_returns_factors,
    const arma::vec& mean_returns
) {

  // Compute sdf coefficients using the FM method
  return SolveSympd(
    covariance_returns_factors.t() * covariance_returns_factors,
    covariance_returns_factors.t()
  ) * mean_returns;

}


////////////////////////////////
///// GKRSDFCoefficientsCpp ////

arma::vec GKRSDFCoefficientsCpp(
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns
) {

  // Solve for `V[R]^{-1} * Cov[R, F]
  const arma::mat var_ret_inv_cov_ret_fac = SolveSympd(
    variance_returns,
    covariance_returns_factors
  );

  // Compute sdf coefficients using the GKR method
  return SolveSympd(
    var_ret_inv_cov_ret_fac.t() * covariance_returns_factors,
    var_ret_inv_cov_ret_fac.t()
  ) * mean_returns;

}

/////////////////////////////////////////////
///// StandardErrorsFMSDFCoefficientsCpp ////

arma::vec StandardErrorsFMSDFCoefficientsCpp(
  const arma::vec& fm_sdf_coefficients,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Compute matrix H used in the standard error calculation
  const arma::mat H = InvSympd(
    covariance_returns_factors.t() * covariance_returns_factors
  );

  // Center returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Term 1: estimation risk for the mean
  const arma::mat term1 = returns_centred * covariance_returns_factors * H;

  // Term 2: covariance
  const arma::mat fac_cen_H = factors_centred * H;
  const arma::mat term2 = fac_cen_H.each_col() % (returns_centred * mean_returns);

  // Term 3 and 4: (covariance * covariance)^{-1}
  const arma::vec cov_ret_fac_mean_ret = covariance_returns_factors.t() * mean_returns;
  const arma::mat ret_cen_cov_ret_fac_H = returns_centred * covariance_returns_factors * H;
  const arma::mat term3 = fac_cen_H.each_col() % (ret_cen_cov_ret_fac_H * cov_ret_fac_mean_ret);
  const arma::mat term4 = ret_cen_cov_ret_fac_H.each_col() % (fac_cen_H * cov_ret_fac_mean_ret);

  // Series
  arma::mat series = term1 + term2 - term3 - term4;

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(series, hac_prewhite) / std::sqrt(
      static_cast<double>(returns.n_rows)
  );

}


// arma::vec StandardErrorsFMSDFCoefficientsCpp(
//   const arma::vec& fm_sdf_coefficients,
//   const arma::mat& returns,
//   const arma::mat& factors,
//   const arma::mat& covariance_returns_factors,
//   const arma::vec& mean_returns,
//   const bool hac_prewhite
// ) {
//
//   // Compute matrix H used in the standard error calculation
//   const arma::mat h_matrix = InvSympd(
//     covariance_returns_factors.t() * covariance_returns_factors
//   );
//
//   // Center returns and factors
//   const arma::mat returns_centred = returns.each_row() - mean_returns.t();
//   const arma::mat factors_centred = factors.each_row() - arma::mean(factors);
//
//   // Term 1: risk due to mean estimation
//   const arma::mat term1 = returns_centred * covariance_returns_factors * h_matrix;
//
//   // Term 2: risk due to covariance estimation
//   const arma::mat psi_centred = term1 - factors_centred;
//   const arma::mat fac_cen_var_fac_inv = SolveSympd(
//     arma::cov(factors), factors_centred.t()
//   ).t();
//   const arma::mat w_t = fac_cen_var_fac_inv * fm_sdf_coefficients;
//   const arma::mat term2 = psi_centred.each_col() % w_t;
//
//   // Term3: misspecification
//   // Compute the errors of the model
//   const arma::vec fm_err = mean_returns -
//     covariance_returns_factors * fm_sdf_coefficients;
//   const arma::mat u_t = returns_centred * fm_err;
//   const arma::mat z_t = fac_cen_var_fac_inv.each_col() % u_t;
//   const arma::mat term3 = z_t * h_matrix;
//
//   // Series
//   arma::mat series = term1 - term2 + term3;
//
//   // Return the HAC standard errors of the estimator
//   return HACStandardErrorsCpp(series, hac_prewhite) / std::sqrt(
//       static_cast<double>(returns.n_rows)
//   );
//
// }


//////////////////////////////////////////////
///// StandardErrorsGKRSDFCoefficientsCpp ////

arma::vec StandardErrorsGKRSDFCoefficientsCpp(
  const arma::vec& gkr_sdf_coefficients,
  const arma::mat& returns,
  const arma::mat& factors,
  const arma::mat& covariance_returns_factors,
  const arma::mat& variance_returns,
  const arma::vec& mean_returns,
  const bool hac_prewhite
) {

  // Compute the inverse of variance_returns and useful matrices that use it
  const arma::mat var_ret_inv = InvSympd(variance_returns);
  const arma::mat var_ret_inv_cov_ret_fac = var_ret_inv * covariance_returns_factors;
  const arma::vec var_ret_inv_mean_ret = var_ret_inv * mean_returns;

  // Compute matrix H used in the standard error calculation
  const arma::mat h_matrix = InvSympd(
    covariance_returns_factors.t() * var_ret_inv_cov_ret_fac
  );

  // Center returns and factors
  const arma::mat returns_centred = returns.each_row() - mean_returns.t();
  const arma::mat factors_centred = factors.each_row() - arma::mean(factors);

  // Term 1
  const arma::mat ARbar = returns_centred * var_ret_inv_cov_ret_fac * h_matrix;
  const arma::vec y =  factors_centred * gkr_sdf_coefficients;
  const arma::mat term1 = ARbar.each_col() % y;

  // Term 2
  const arma::vec err = mean_returns - covariance_returns_factors * gkr_sdf_coefficients;
  const arma::vec u = returns_centred * SolveSympd(variance_returns, err);
  const arma::mat HFbarMinusARbar = factors_centred * h_matrix - ARbar;
  const arma::mat term2 = HFbarMinusARbar.each_col() % u;

  // all terms
  arma::mat series = term1 + term2 + arma::repmat(
    gkr_sdf_coefficients.t(), returns.n_rows, 1
  );

  // Return the HAC standard errors of the estimator
  return HACStandardErrorsCpp(series, hac_prewhite) / std::sqrt(
      static_cast<double>(returns.n_rows)
  );

}


// arma::vec StandardErrorsGKRSDFCoefficientsCpp(
//   const arma::vec& gkr_sdf_coefficients,
//   const arma::mat& returns,
//   const arma::mat& factors,
//   const arma::mat& covariance_returns_factors,
//   const arma::mat& variance_returns,
//   const arma::vec& mean_returns,
//   const bool hac_prewhite
// ) {
//
//   // Compute the inverse of variance_returns and useful matrices that use it
//   const arma::mat var_ret_inv = InvSympd(variance_returns);
//   const arma::mat var_ret_inv_cov_ret_fac = var_ret_inv * covariance_returns_factors;
//   const arma::vec var_ret_inv_mean_ret = var_ret_inv * mean_returns;
//
//   // Compute matrix H used in the standard error calculation
//   const arma::mat h_matrix = InvSympd(
//     covariance_returns_factors.t() * var_ret_inv_cov_ret_fac
//   );
//
//   // Center returns and factors
//   const arma::mat returns_centred = returns.each_row() - mean_returns.t();
//   const arma::mat factors_centred = factors.each_row() - arma::mean(factors);
//
//   // Term 1: risk due to mean estimation
//   const arma::mat term1 = returns_centred * var_ret_inv_cov_ret_fac * h_matrix;
//
//   // Term 2: risk due to covariance estimation
//   const arma::mat psi_centred = term1 - factors_centred;
//   const arma::mat fac_cen_var_fac_inv = SolveSympd(
//     arma::cov(factors), factors_centred.t()
//   ).t();
//   const arma::mat w_t = fac_cen_var_fac_inv * gkr_sdf_coefficients;
//   const arma::mat term2 = psi_centred.each_col() % w_t;
//
//   // Term3: misspecification1
//   // Compute the errors of the model
//   const arma::vec gkr_err = mean_returns -
//     covariance_returns_factors * gkr_sdf_coefficients;
//   const arma::mat u_t = returns_centred * SolveSympd(
//     variance_returns, gkr_err
//   );
//   const arma::mat z_t = fac_cen_var_fac_inv.each_col() % u_t;
//   const arma::mat term3 = z_t * h_matrix;
//
//   // Term4: misspecification2
//   const arma::mat term4 = term1.each_col() % u_t;
//
//   // Series
//   arma::mat series = term1 - term2 + term3 - term4;
//
//   // Return the HAC standard errors of the estimator
//   return HACStandardErrorsCpp(series, hac_prewhite) / std::sqrt(
//       static_cast<double>(returns.n_rows)
//   );
//
// }

// arma::vec StandardErrorsGKRSDFCoefficientsCpp1(
//   const arma::vec& gkr_sdf_coefficients,
//   const arma::mat& returns,
//   const arma::mat& factors,
//   const arma::mat& covariance_returns_factors,
//   const arma::mat& variance_returns,
//   const arma::vec& mean_returns,
//   const bool hac_prewhite
// ) {
//
//   // Compute the inverse of variance_returns and useful matrices that use it
//   const arma::mat var_ret_inv = InvSympd(variance_returns);
//   const arma::mat var_ret_inv_cov_ret_fac = var_ret_inv * covariance_returns_factors;
//   const arma::vec var_ret_inv_mean_ret = var_ret_inv * mean_returns;
//
//   // Compute matrix H used in the standard error calculation
//   const arma::mat h_matrix = InvSympd(
//     covariance_returns_factors.t() * var_ret_inv_cov_ret_fac
//   );
//
//   // Center returns and factors
//   const arma::mat returns_centred = returns.each_row() - mean_returns.t();
//   const arma::mat factors_centred = factors.each_row() - arma::mean(factors);
//
//   // Compute intermediate terms for standard error calculation
//   const arma::mat fac_cen_gkr_sdf_coeff = factors_centred * gkr_sdf_coefficients;
//   const arma::mat ret_cen_var_ret_inv = returns_centred * var_ret_inv;
//   const arma::mat ret_cen_fac_cen_gkr_sdf_coeff = returns_centred.each_col() % fac_cen_gkr_sdf_coeff;
//   const arma::mat gkr_err = ret_cen_fac_cen_gkr_sdf_coeff.each_row() - mean_returns.t();
//   const arma::mat mimicking_err_h =
//     (factors_centred - ret_cen_var_ret_inv * covariance_returns_factors) *
//     h_matrix;
//   const arma::vec u = arma::sum(gkr_err % ret_cen_var_ret_inv, 1);
//
//   // Compute terms used in the standard error calculation
//   const arma::mat term1 = gkr_err * var_ret_inv_cov_ret_fac * h_matrix;
//   const arma::mat term2 = mimicking_err_h.each_col() % u;
//
//   arma::mat series = term1 + term2;
//
//   // Return the HAC standard errors of the estimator
//   return HACStandardErrorsCpp(series, hac_prewhite) / std::sqrt(
//       static_cast<double>(returns.n_rows)
//   );
//
// }
