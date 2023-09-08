// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// OptimalAdaptiveFRPGCVCpp
Rcpp::List OptimalAdaptiveFRPGCVCpp(const arma::mat& returns, const arma::mat& factors, const arma::mat& covariance_factors_returns, const arma::mat& variance_returns, const arma::vec& mean_returns, const arma::vec& penalty_parameters, const char weighting_type, const bool one_stddev_rule, const bool gcv_vr_weighting, const bool gcv_scaling_n_assets, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_OptimalAdaptiveFRPGCVCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP covariance_factors_returnsSEXP, SEXP variance_returnsSEXP, SEXP mean_returnsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP gcv_vr_weightingSEXP, SEXP gcv_scaling_n_assetsSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance_factors_returns(covariance_factors_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type variance_returns(variance_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean_returns(mean_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const bool >::type gcv_vr_weighting(gcv_vr_weightingSEXP);
    Rcpp::traits::input_parameter< const bool >::type gcv_scaling_n_assets(gcv_scaling_n_assetsSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(OptimalAdaptiveFRPGCVCpp(returns, factors, covariance_factors_returns, variance_returns, mean_returns, penalty_parameters, weighting_type, one_stddev_rule, gcv_vr_weighting, gcv_scaling_n_assets, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// OptimalAdaptiveFRPCVCpp
Rcpp::List OptimalAdaptiveFRPCVCpp(const arma::mat& returns, const arma::mat& factors, const arma::mat& covariance_factors_returns, const arma::mat& variance_returns, const arma::vec& mean_returns, const arma::vec& penalty_parameters, const char weighting_type, const bool one_stddev_rule, const unsigned int n_folds, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_OptimalAdaptiveFRPCVCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP covariance_factors_returnsSEXP, SEXP variance_returnsSEXP, SEXP mean_returnsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP n_foldsSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance_factors_returns(covariance_factors_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type variance_returns(variance_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean_returns(mean_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_folds(n_foldsSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(OptimalAdaptiveFRPCVCpp(returns, factors, covariance_factors_returns, variance_returns, mean_returns, penalty_parameters, weighting_type, one_stddev_rule, n_folds, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// OptimalAdaptiveFRPRVCpp
Rcpp::List OptimalAdaptiveFRPRVCpp(const arma::mat& returns, const arma::mat& factors, const arma::mat& covariance_factors_returns, const arma::mat& variance_returns, const arma::vec& mean_returns, const arma::vec& penalty_parameters, const char weighting_type, const bool one_stddev_rule, const unsigned int n_train_observations, const unsigned int n_test_observations, const unsigned int roll_shift, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_OptimalAdaptiveFRPRVCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP covariance_factors_returnsSEXP, SEXP variance_returnsSEXP, SEXP mean_returnsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP n_train_observationsSEXP, SEXP n_test_observationsSEXP, SEXP roll_shiftSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance_factors_returns(covariance_factors_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type variance_returns(variance_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean_returns(mean_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_train_observations(n_train_observationsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_test_observations(n_test_observationsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type roll_shift(roll_shiftSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(OptimalAdaptiveFRPRVCpp(returns, factors, covariance_factors_returns, variance_returns, mean_returns, penalty_parameters, weighting_type, one_stddev_rule, n_train_observations, n_test_observations, roll_shift, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// OptimalAdaptiveIFRPGCVCpp
Rcpp::List OptimalAdaptiveIFRPGCVCpp(const arma::mat& returns, const arma::mat& factors, const arma::mat& covariance_factors_returns, const arma::mat& variance_returns, const arma::vec& mean_returns, const arma::vec& penalty_parameters, const char weighting_type, const bool one_stddev_rule, const bool gcv_vr_weighting, const bool gcv_scaling_n_assets, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_OptimalAdaptiveIFRPGCVCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP covariance_factors_returnsSEXP, SEXP variance_returnsSEXP, SEXP mean_returnsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP gcv_vr_weightingSEXP, SEXP gcv_scaling_n_assetsSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance_factors_returns(covariance_factors_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type variance_returns(variance_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean_returns(mean_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const bool >::type gcv_vr_weighting(gcv_vr_weightingSEXP);
    Rcpp::traits::input_parameter< const bool >::type gcv_scaling_n_assets(gcv_scaling_n_assetsSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(OptimalAdaptiveIFRPGCVCpp(returns, factors, covariance_factors_returns, variance_returns, mean_returns, penalty_parameters, weighting_type, one_stddev_rule, gcv_vr_weighting, gcv_scaling_n_assets, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// OptimalAdaptiveIFRPCVCpp
Rcpp::List OptimalAdaptiveIFRPCVCpp(const arma::mat& returns, const arma::mat& factors, const arma::mat& covariance_factors_returns, const arma::mat& variance_returns, const arma::vec& mean_returns, const arma::vec& penalty_parameters, const char weighting_type, const bool one_stddev_rule, const unsigned int n_folds, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_OptimalAdaptiveIFRPCVCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP covariance_factors_returnsSEXP, SEXP variance_returnsSEXP, SEXP mean_returnsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP n_foldsSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance_factors_returns(covariance_factors_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type variance_returns(variance_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean_returns(mean_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_folds(n_foldsSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(OptimalAdaptiveIFRPCVCpp(returns, factors, covariance_factors_returns, variance_returns, mean_returns, penalty_parameters, weighting_type, one_stddev_rule, n_folds, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// OptimalAdaptiveIFRPRVCpp
Rcpp::List OptimalAdaptiveIFRPRVCpp(const arma::mat& returns, const arma::mat& factors, const arma::mat& covariance_factors_returns, const arma::mat& variance_returns, const arma::vec& mean_returns, const arma::vec& penalty_parameters, const char weighting_type, const bool one_stddev_rule, const unsigned int n_train_observations, const unsigned int n_test_observations, const unsigned int roll_shift, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_OptimalAdaptiveIFRPRVCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP covariance_factors_returnsSEXP, SEXP variance_returnsSEXP, SEXP mean_returnsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP n_train_observationsSEXP, SEXP n_test_observationsSEXP, SEXP roll_shiftSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type covariance_factors_returns(covariance_factors_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type variance_returns(variance_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type mean_returns(mean_returnsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_train_observations(n_train_observationsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_test_observations(n_test_observationsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type roll_shift(roll_shiftSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(OptimalAdaptiveIFRPRVCpp(returns, factors, covariance_factors_returns, variance_returns, mean_returns, penalty_parameters, weighting_type, one_stddev_rule, n_train_observations, n_test_observations, roll_shift, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// AdaptiveIFRPCpp
arma::mat AdaptiveIFRPCpp(const arma::vec& ifrp, const arma::vec& weights, const arma::vec& penalty_parameter);
RcppExport SEXP _intrinsicFRP_AdaptiveIFRPCpp(SEXP ifrpSEXP, SEXP weightsSEXP, SEXP penalty_parameterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type ifrp(ifrpSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameter(penalty_parameterSEXP);
    rcpp_result_gen = Rcpp::wrap(AdaptiveIFRPCpp(ifrp, weights, penalty_parameter));
    return rcpp_result_gen;
END_RCPP
}
// AdaptiveWeightsCpp
arma::vec AdaptiveWeightsCpp(const arma::mat& returns, const arma::mat& factors, const char type);
RcppExport SEXP _intrinsicFRP_AdaptiveWeightsCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const char >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(AdaptiveWeightsCpp(returns, factors, type));
    return rcpp_result_gen;
END_RCPP
}
// FRPCpp
Rcpp::List FRPCpp(const arma::mat& returns, const arma::mat& factors, const bool misspecification_robust, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_FRPCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP misspecification_robustSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type misspecification_robust(misspecification_robustSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(FRPCpp(returns, factors, misspecification_robust, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}
// IFRPCpp
Rcpp::List IFRPCpp(const arma::mat& returns, const arma::mat& factors, const bool include_standard_errors);
RcppExport SEXP _intrinsicFRP_IFRPCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP include_standard_errorsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    rcpp_result_gen = Rcpp::wrap(IFRPCpp(returns, factors, include_standard_errors));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_intrinsicFRP_OptimalAdaptiveFRPGCVCpp", (DL_FUNC) &_intrinsicFRP_OptimalAdaptiveFRPGCVCpp, 11},
    {"_intrinsicFRP_OptimalAdaptiveFRPCVCpp", (DL_FUNC) &_intrinsicFRP_OptimalAdaptiveFRPCVCpp, 10},
    {"_intrinsicFRP_OptimalAdaptiveFRPRVCpp", (DL_FUNC) &_intrinsicFRP_OptimalAdaptiveFRPRVCpp, 12},
    {"_intrinsicFRP_OptimalAdaptiveIFRPGCVCpp", (DL_FUNC) &_intrinsicFRP_OptimalAdaptiveIFRPGCVCpp, 11},
    {"_intrinsicFRP_OptimalAdaptiveIFRPCVCpp", (DL_FUNC) &_intrinsicFRP_OptimalAdaptiveIFRPCVCpp, 10},
    {"_intrinsicFRP_OptimalAdaptiveIFRPRVCpp", (DL_FUNC) &_intrinsicFRP_OptimalAdaptiveIFRPRVCpp, 12},
    {"_intrinsicFRP_AdaptiveIFRPCpp", (DL_FUNC) &_intrinsicFRP_AdaptiveIFRPCpp, 3},
    {"_intrinsicFRP_AdaptiveWeightsCpp", (DL_FUNC) &_intrinsicFRP_AdaptiveWeightsCpp, 3},
    {"_intrinsicFRP_FRPCpp", (DL_FUNC) &_intrinsicFRP_FRPCpp, 4},
    {"_intrinsicFRP_IFRPCpp", (DL_FUNC) &_intrinsicFRP_IFRPCpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_intrinsicFRP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
