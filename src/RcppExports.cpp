// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// FGXThreePassCovarianceCpp
arma::mat FGXThreePassCovarianceCpp(const arma::mat& returns, const arma::mat& control_factors, const arma::mat& new_factors, const arma::vec& sdf_coefficients, const arma::uvec& idx_selected);
RcppExport SEXP _intrinsicFRP_FGXThreePassCovarianceCpp(SEXP returnsSEXP, SEXP control_factorsSEXP, SEXP new_factorsSEXP, SEXP sdf_coefficientsSEXP, SEXP idx_selectedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type control_factors(control_factorsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type new_factors(new_factorsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type sdf_coefficients(sdf_coefficientsSEXP);
    Rcpp::traits::input_parameter< const arma::uvec& >::type idx_selected(idx_selectedSEXP);
    rcpp_result_gen = Rcpp::wrap(FGXThreePassCovarianceCpp(returns, control_factors, new_factors, sdf_coefficients, idx_selected));
    return rcpp_result_gen;
END_RCPP
}
// FRPCpp
Rcpp::List FRPCpp(const arma::mat& returns, const arma::mat& factors, const bool misspecification_robust, const bool include_standard_errors, const bool hac_prewhite, const double target_level_gkr2014_screening);
RcppExport SEXP _intrinsicFRP_FRPCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP misspecification_robustSEXP, SEXP include_standard_errorsSEXP, SEXP hac_prewhiteSEXP, SEXP target_level_gkr2014_screeningSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type misspecification_robust(misspecification_robustSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type hac_prewhite(hac_prewhiteSEXP);
    Rcpp::traits::input_parameter< const double >::type target_level_gkr2014_screening(target_level_gkr2014_screeningSEXP);
    rcpp_result_gen = Rcpp::wrap(FRPCpp(returns, factors, misspecification_robust, include_standard_errors, hac_prewhite, target_level_gkr2014_screening));
    return rcpp_result_gen;
END_RCPP
}
// GKRFactorScreeningCpp
Rcpp::List GKRFactorScreeningCpp(const arma::mat& returns, const arma::mat& factors, const double target_level, const bool hac_prewhite);
RcppExport SEXP _intrinsicFRP_GKRFactorScreeningCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP target_levelSEXP, SEXP hac_prewhiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const double >::type target_level(target_levelSEXP);
    Rcpp::traits::input_parameter< const bool >::type hac_prewhite(hac_prewhiteSEXP);
    rcpp_result_gen = Rcpp::wrap(GKRFactorScreeningCpp(returns, factors, target_level, hac_prewhite));
    return rcpp_result_gen;
END_RCPP
}
// HACCovarianceMatrixCpp
arma::mat HACCovarianceMatrixCpp(arma::mat& series, const bool prewhite);
RcppExport SEXP _intrinsicFRP_HACCovarianceMatrixCpp(SEXP seriesSEXP, SEXP prewhiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat& >::type series(seriesSEXP);
    Rcpp::traits::input_parameter< const bool >::type prewhite(prewhiteSEXP);
    rcpp_result_gen = Rcpp::wrap(HACCovarianceMatrixCpp(series, prewhite));
    return rcpp_result_gen;
END_RCPP
}
// HACVarianceCpp
double HACVarianceCpp(arma::vec& series, const bool prewhite);
RcppExport SEXP _intrinsicFRP_HACVarianceCpp(SEXP seriesSEXP, SEXP prewhiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec& >::type series(seriesSEXP);
    Rcpp::traits::input_parameter< const bool >::type prewhite(prewhiteSEXP);
    rcpp_result_gen = Rcpp::wrap(HACVarianceCpp(series, prewhite));
    return rcpp_result_gen;
END_RCPP
}
// HJMisspecificationDistanceCpp
Rcpp::List HJMisspecificationDistanceCpp(const arma::mat& returns, const arma::mat& factors, const double ci_coverage, const bool hac_prewhite);
RcppExport SEXP _intrinsicFRP_HJMisspecificationDistanceCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP ci_coverageSEXP, SEXP hac_prewhiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const double >::type ci_coverage(ci_coverageSEXP);
    Rcpp::traits::input_parameter< const bool >::type hac_prewhite(hac_prewhiteSEXP);
    rcpp_result_gen = Rcpp::wrap(HJMisspecificationDistanceCpp(returns, factors, ci_coverage, hac_prewhite));
    return rcpp_result_gen;
END_RCPP
}
// ChenFang2019BetaRankTestCpp
Rcpp::List ChenFang2019BetaRankTestCpp(const arma::mat& returns, const arma::mat& factors, const unsigned int n_bootstrap, const double target_level_kp2006_rank_test);
RcppExport SEXP _intrinsicFRP_ChenFang2019BetaRankTestCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP n_bootstrapSEXP, SEXP target_level_kp2006_rank_testSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_bootstrap(n_bootstrapSEXP);
    Rcpp::traits::input_parameter< const double >::type target_level_kp2006_rank_test(target_level_kp2006_rank_testSEXP);
    rcpp_result_gen = Rcpp::wrap(ChenFang2019BetaRankTestCpp(returns, factors, n_bootstrap, target_level_kp2006_rank_test));
    return rcpp_result_gen;
END_RCPP
}
// IterativeKleibergenPaap2006BetaRankTestCpp
Rcpp::List IterativeKleibergenPaap2006BetaRankTestCpp(const arma::mat& returns, const arma::mat& factors, const double target_level);
RcppExport SEXP _intrinsicFRP_IterativeKleibergenPaap2006BetaRankTestCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP target_levelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const double >::type target_level(target_levelSEXP);
    rcpp_result_gen = Rcpp::wrap(IterativeKleibergenPaap2006BetaRankTestCpp(returns, factors, target_level));
    return rcpp_result_gen;
END_RCPP
}
// OracleTFRPCpp
Rcpp::List OracleTFRPCpp(const arma::mat& returns, const arma::mat& factors, const arma::vec& penalty_parameters, const char weighting_type, const char tuning_type, const bool one_stddev_rule, const bool gcv_scaling_n_assets, const bool gcv_identification_check, const double target_level_kp2006_rank_test, const unsigned int n_folds, const unsigned int n_train_observations, const unsigned int n_test_observations, const unsigned int roll_shift, const bool relaxed, const bool include_standard_errors, const bool hac_prewhite);
RcppExport SEXP _intrinsicFRP_OracleTFRPCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP penalty_parametersSEXP, SEXP weighting_typeSEXP, SEXP tuning_typeSEXP, SEXP one_stddev_ruleSEXP, SEXP gcv_scaling_n_assetsSEXP, SEXP gcv_identification_checkSEXP, SEXP target_level_kp2006_rank_testSEXP, SEXP n_foldsSEXP, SEXP n_train_observationsSEXP, SEXP n_test_observationsSEXP, SEXP roll_shiftSEXP, SEXP relaxedSEXP, SEXP include_standard_errorsSEXP, SEXP hac_prewhiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type penalty_parameters(penalty_parametersSEXP);
    Rcpp::traits::input_parameter< const char >::type weighting_type(weighting_typeSEXP);
    Rcpp::traits::input_parameter< const char >::type tuning_type(tuning_typeSEXP);
    Rcpp::traits::input_parameter< const bool >::type one_stddev_rule(one_stddev_ruleSEXP);
    Rcpp::traits::input_parameter< const bool >::type gcv_scaling_n_assets(gcv_scaling_n_assetsSEXP);
    Rcpp::traits::input_parameter< const bool >::type gcv_identification_check(gcv_identification_checkSEXP);
    Rcpp::traits::input_parameter< const double >::type target_level_kp2006_rank_test(target_level_kp2006_rank_testSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_folds(n_foldsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_train_observations(n_train_observationsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type n_test_observations(n_test_observationsSEXP);
    Rcpp::traits::input_parameter< const unsigned int >::type roll_shift(roll_shiftSEXP);
    Rcpp::traits::input_parameter< const bool >::type relaxed(relaxedSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type hac_prewhite(hac_prewhiteSEXP);
    rcpp_result_gen = Rcpp::wrap(OracleTFRPCpp(returns, factors, penalty_parameters, weighting_type, tuning_type, one_stddev_rule, gcv_scaling_n_assets, gcv_identification_check, target_level_kp2006_rank_test, n_folds, n_train_observations, n_test_observations, roll_shift, relaxed, include_standard_errors, hac_prewhite));
    return rcpp_result_gen;
END_RCPP
}
// TFRPCpp
Rcpp::List TFRPCpp(const arma::mat& returns, const arma::mat& factors, const bool include_standard_errors, const bool hac_prewhite);
RcppExport SEXP _intrinsicFRP_TFRPCpp(SEXP returnsSEXP, SEXP factorsSEXP, SEXP include_standard_errorsSEXP, SEXP hac_prewhiteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type returns(returnsSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type factors(factorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type include_standard_errors(include_standard_errorsSEXP);
    Rcpp::traits::input_parameter< const bool >::type hac_prewhite(hac_prewhiteSEXP);
    rcpp_result_gen = Rcpp::wrap(TFRPCpp(returns, factors, include_standard_errors, hac_prewhite));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_intrinsicFRP_FGXThreePassCovarianceCpp", (DL_FUNC) &_intrinsicFRP_FGXThreePassCovarianceCpp, 5},
    {"_intrinsicFRP_FRPCpp", (DL_FUNC) &_intrinsicFRP_FRPCpp, 6},
    {"_intrinsicFRP_GKRFactorScreeningCpp", (DL_FUNC) &_intrinsicFRP_GKRFactorScreeningCpp, 4},
    {"_intrinsicFRP_HACCovarianceMatrixCpp", (DL_FUNC) &_intrinsicFRP_HACCovarianceMatrixCpp, 2},
    {"_intrinsicFRP_HACVarianceCpp", (DL_FUNC) &_intrinsicFRP_HACVarianceCpp, 2},
    {"_intrinsicFRP_HJMisspecificationDistanceCpp", (DL_FUNC) &_intrinsicFRP_HJMisspecificationDistanceCpp, 4},
    {"_intrinsicFRP_ChenFang2019BetaRankTestCpp", (DL_FUNC) &_intrinsicFRP_ChenFang2019BetaRankTestCpp, 4},
    {"_intrinsicFRP_IterativeKleibergenPaap2006BetaRankTestCpp", (DL_FUNC) &_intrinsicFRP_IterativeKleibergenPaap2006BetaRankTestCpp, 3},
    {"_intrinsicFRP_OracleTFRPCpp", (DL_FUNC) &_intrinsicFRP_OracleTFRPCpp, 16},
    {"_intrinsicFRP_TFRPCpp", (DL_FUNC) &_intrinsicFRP_TFRPCpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_intrinsicFRP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
