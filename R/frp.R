# Author: Alberto Quaini

###################################
######  Factor Risk Premia ########
###################################

#' @title Factor risk premia.
#'
#' @name FRP
#' @description Computes the Fama-MachBeth (1973) <doi:10.1086/260061> factor
#' risk premia:
#' `FMFRP = (beta' * beta)^{-1} * beta' * E[R]` where
#' `beta = Cov[R, F] * V[F]^{-1}`
#' or the misspecification-robust factor risk premia of Kan-Robotti-Shanken (2013)
#' <doi:10.1111/jofi.12035>:
#' `KRSFRP = (beta' * V[R]^{-1} * beta)^{-1} * beta' * V[R]^{-1} * E[R]`,
#' from data on factors `F` and test
#' asset excess returns `R`.
#' These notions of factor risk premia are by construction the negative
#' covariance of factors `F` with candidate SDF
#' `M = 1 - d' * (F - E[F])`,
#' where SDF coefficients `d` are obtained by minimizing pricing errors:
#' `argmin_{d} (E[R] - Cov[R,F] * d)' * (E[R] - Cov[R,F] * d)`
#' and
#' `argmin_{d} (E[R] - Cov[R,F] * d)' * V[R]^{-1} * (E[R] - Cov[R,F] * d)`,
#' respectively.
#' Optionally computes the corresponding
#' heteroskedasticity and autocorrelation robust standard errors (accounting
#' for a potential model misspecification) using the
#' Newey-West (1994) <doi:10.2307/2297912> plug-in procedure to select the
#' number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
#' For the standard error computations, the function allows to internally
#' pre-whiten the series by fitting a VAR(1),
#' i.e., a vector autoregressive model of order 1.
#' All the details can be found in Kan-Robotti-Shanken (2013)
#' <doi:10.1111/jofi.12035>.
#'
#' @param returns A `n_observations x n_returns`-dimensional matrix of test asset
#' excess returns.
#' @param factors A `n_observations x n_factors`-dimensional matrix of factors.
#' @param misspecification_robust A boolean: `TRUE` for the
#' "misspecification-robust" Kan-Robotti-Shanken (2013) GLS approach using the
#' inverse covariance matrix of returns; `FALSE` for standard Fama-MacBeth
#' risk premia. Default is `TRUE`.
#' @param include_standard_errors A boolean: `TRUE` if you want to compute the
#' factor risk premia HAC standard errors; `FALSE` otherwise.
#' Default is `FALSE`.
#' @param hac_prewhite A boolean indicating if the series needs pre-whitening by
#' fitting an AR(1) in the internal heteroskedasticity and autocorrelation
#' robust covariance (HAC) estimation. Default is `false`.
#' @param target_level_gkr2014_screening A number indicating the target level of
#' the tests underlying the factor screening procedure in Gospodinov-Kan-Robotti
#' (2014). If it is zero, then no factor screening procedure is
#' implemented. Otherwise, it implements an iterative screening procedure
#' based on the sequential removal of factors associated with the smallest insignificant
#' t-test of a nonzero SDF coefficient. The threshold for the absolute t-test is
#' `target_level_gkr2014_screening / n_factors`, where n_factors indicate the
#' number of factors in the model at the current iteration. Default is `0.`, i.e.,
#' no factor screening.
#' @param check_arguments A boolean: `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A list containing `n_factors`-dimensional vector of factor
#' risk premia in `"risk_premia"`; if `include_standard_errors = TRUE`, then
#' it further includes `n_factors`-dimensional vector of factor risk
#' premia standard errors in `"standard_errors"`;
#' if `target_level_gkr2014_screening >= 0`, it further includes the indices of
#' the selected factors in `selected_factor_indices`.
#'
#' @examples
#' # import package data on 6 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' # compute KRS factor risk premia and their standard errors
#' frp = FRP(returns, factors, include_standard_errors = TRUE)
#'
#' @export
FRP = function(
  returns,
  factors,
  misspecification_robust = TRUE,
  include_standard_errors = FALSE,
  hac_prewhite = FALSE,
  target_level_gkr2014_screening = 0.,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {

    CheckData(returns, factors)
    stopifnot("`misspecification_robust` must be boolean" = is.logical(misspecification_robust))
    stopifnot("`include_standard_errors` must be boolean" = is.logical(include_standard_errors))
    stopifnot("`hac_prewhite` must be boolean" = is.logical(hac_prewhite))
    stopifnot("`target_level_gkr2014_screening` must be numeric" = is.numeric(target_level_gkr2014_screening))
    stopifnot("`target_level_gkr2014_screening` must be between 0 and 1" = (target_level_gkr2014_screening >= 0.) & (target_level_gkr2014_screening <= 1.))

  }

  # compute the FRP estimate and, eventually, their standard errors
  return(.Call(`_intrinsicFRP_FRPCpp`,
    returns,
    factors,
    misspecification_robust,
    include_standard_errors,
    hac_prewhite,
    target_level_gkr2014_screening
  ))

}

#' @title Compute Factor Risk Premia using Giglio and Xiu (2021) Method
#'
#' @description Computes the factor risk premia using the 3-step procedure of
#' Giglio and Xiu (2021) <doi:10.1086/714090>. The procedure consists in
#' 1) extracting p PCA from returns, where p is a consistent estimator of the
#' number of strong latent factors, 2) compute a cross-sectional regression of
#' average returns on the estimated betas of the p latent factors, and 3) compute
#' a time-series regression of observed factors on the p latent factors to obtain
#' their mimicking portfolio risk premia.
#'
#' @param returns A `n_observations x n_returns` matrix of asset returns.
#' @param factors A `n_observations x n_factors` matrix of risk factors.
#' @param which_n_pca Method to determine the number of PCA components:
#' `0` for Giglio Xiu (2021) <doi:10.1086/714090>,
#' `-1` for Ahn Horenstein (2013) <doi:10.3982/ECTA8968>,
#' or any positive integer for a specific number of PCs.
#' @param check_arguments A boolean: `TRUE` for internal check of all function
#' arguments; `FALSE` otherwise. Default is `TRUE`.
#'
#' @return A list containing:
#' \describe{
#'   \item{risk_premia}{A matrix of factor risk premia.}
#'   \item{n_pca}{The number of principal components used.}
#' }
#' @examples
#' \dontrun{
#' returns <- matrix(rnorm(200), nrow=20, ncol=10)
#' factors <- matrix(rnorm(60), nrow=20, ncol=3)
#' which_n_pca <- 0
#' result <- GiglioXiu2021RiskPremia(returns, factors, which_n_pca)
#' print(result)
#' }
#' @export
GiglioXiu2021RiskPremia = function(
  returns,
  factors,
  which_n_pca = 0,
  check_arguments = TRUE
) {

  # check function arguments
  if (check_arguments) {

    stopifnot("`returns` must contain numeric values" = is.numeric(returns))
    stopifnot("`factors` must contain numeric values" = is.numeric(factors))
    stopifnot("`returns` and `factors` must have the same number of rows" = nrow(returns) == nrow(factors))
    stopifnot("`returns` must not contain missing values (NA/NaN)" = !anyNA(returns))
    stopifnot("`factors` must not contain missing values (NA/NaN)" = !anyNA(factors))

    stopifnot("`which_n_pca` must be numeric" = is.numeric(which_n_pca))

  }

  # Ensure the function is defined in the C++ file and compiled
  return(.Call(`_intrinsicFRP_GiglioXiu2021RiskPremiaCpp`,
    returns,
    factors,
    which_n_pca
  ))

}
