% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oracle_tfrp.R
\name{OracleTFRP}
\alias{OracleTFRP}
\title{Oracle tradable factor risk premia.}
\usage{
OracleTFRP(
  returns,
  factors,
  penalty_parameters,
  weighting_type = "c",
  tuning_type = "g",
  one_stddev_rule = TRUE,
  gcv_scaling_n_assets = FALSE,
  gcv_identification_check = FALSE,
  target_level_kp2006_rank_test = 0.05,
  n_folds = 5,
  n_train_observations = 120,
  n_test_observations = 12,
  roll_shift = 12,
  relaxed = FALSE,
  include_standard_errors = FALSE,
  hac_prewhite = FALSE,
  plot_score = TRUE,
  check_arguments = TRUE
)
}
\arguments{
\item{returns}{A \verb{n_observations x n_returns}-dimensional matrix of test asset
excess returns.}

\item{factors}{A \verb{n_observations x n_factors}-dimensional matrix of factors.}

\item{penalty_parameters}{A \code{n_parameters}-dimensional vector of penalty
parameter values from smallest to largest.}

\item{weighting_type}{A character specifying the type of adaptive weights:
based on the correlation between factors and returns \code{'c'}; based on the
regression coefficients of returns on factors \code{'b'}; based on the first-step
tradable risk premia estimator \code{'a'}; otherwise a vector of ones (any other
character). Default is \code{'c'}.}

\item{tuning_type}{A character indicating the parameter tuning type:
\code{'g'} for generalized cross validation; \code{'c'} for K-fold cross validation;
\code{'r'} for rolling validation. Default is \code{'g'}.}

\item{one_stddev_rule}{A boolean: \code{TRUE} for picking the most parsimonious model
whose score is not higher than one standard error above the score of the
best model; \code{FALSE} for picking the best model. Default is \code{TRUE}.}

\item{gcv_scaling_n_assets}{(only relevant for \code{tuning_type ='g'})
A boolean: \code{TRUE} for sqrt(n_assets) scaling (\code{sqrt(n_assets) / n_observations});
\code{FALSE} otherwise (\code{1 / n_observations}). Default is \code{FALSE}.}

\item{gcv_identification_check}{(only relevant for \code{tuning_type ='g'})
A boolean: \code{TRUE} for a loose check for model
identification; \code{FALSE} otherwise. Default is \code{FALSE}.}

\item{target_level_kp2006_rank_test}{(only relevant for \code{tuning_type ='g'}
and if \code{gcv_identification_check} is
\code{TRUE}) A number indicating the level of the Kleibergen Paap 2006 rank test.
If it is strictly grater than zero, then the iterative Kleibergen Paap 2006 rank
test at \code{level = target_level_kp2006_rank_test / n_factors}
(where division by the number of factors is performed as a Bonferroni correction for
multiple testing) is used to compute an initial estimator
of the rank of the factor loadings in the Chen Fang 2019 rank test.
Otherwise, the initial rank estimator is taken to be the number of singular
values above \code{n_observations^(-1/4)}.
Default is \code{0.05}.}

\item{n_folds}{(only relevant for \code{tuning_type ='c'}) An integer indicating
the number of k-fold for cross validation. Default is \code{5}.}

\item{n_train_observations}{(only relevant for \code{tuning_type ='r'}) The number of
observations in the rolling training set. Default is \code{120}.}

\item{n_test_observations}{(only relevant for \code{tuning_type ='r'}) The number of
observations in the test set. Default is \code{12}.}

\item{roll_shift}{(only relevant for \code{tuning_type ='r'}) The number of observation
shift when moving from the rolling window to the next one. Default is \code{12}.}

\item{relaxed}{A boolean: \code{TRUE} if you want to compute a post-selection
unpenalized tradable factor risk premia to remove the bias due to shrinkage;
FALSE\verb{otherwise. Default is}FALSE`.}

\item{include_standard_errors}{A boolean \code{TRUE} if you want to compute the
adaptive tradable factor risk premia HAC standard errors; \code{FALSE} otherwise.
Default is \code{FALSE}.}

\item{hac_prewhite}{A boolean indicating if the series needs prewhitening by
fitting an AR(1) in the internal heteroskedasticity and autocorrelation
robust covariance (HAC) estimation. Default is \code{false}.}

\item{plot_score}{A boolean: \code{TRUE} for plotting the model score; \code{FALSE} otherwise.
Default is \code{TRUE}.}

\item{check_arguments}{A boolean \code{TRUE} for internal check of all function
arguments; \code{FALSE} otherwise. Default is \code{TRUE}.}
}
\value{
A list containing the \code{n_factors}-dimensional vector of adaptive
tradable factor risk premia in \code{"risk_premia"}; the optimal penalty
parameter value in \code{"penalty_parameter"}; the model score for each penalty
parameter value in \code{"model_score"};  if \code{include_standard_errors = TRUE}, then
it further includes \code{n_factors}-dimensional vector of tradable factor risk
premia standard errors in \code{"standard_errors"}.
}
\description{
Computes Oracle tradable factor risk premia of
Quaini-Trojani-Yuan (2023) \url{doi:10.2139/ssrn.4574683} from data on
\code{K} factors \verb{F = [F_1,...,F_K]'} and test asset excess returns \code{R}:
\verb{OTFRP = argmin_x ||TFRP - x||_2^2 + tau * sum_\{k=1\}^K w_k * |x_k|},
where \code{TFRP} is the tradable factor risk premia estimator, \code{tau > 0} is a
penalty parameter, and the Oracle weights are given by
\verb{w_k = 1 / ||corr[F_k, R]||_2^2}.
This estimator is called "Oracle" in the sense that the probability that
the index set of its nonzero estimated risk premia equals the index set of
the true strong factors tends to 1 (Oracle selection), and that on the strong
factors, the estimator reaches the optimal asymptotic Normal distribution.
Here, strong factors are those that have a nonzero population marginal
correlation with asset excess returns.
Tuning of the penalty parameter \code{tau} is performed via Generalized Cross
Validation (GCV), Cross Validation (CV) or Rolling Validation (RV).
GCV tunes parameter \code{tau} by minimizing the criterium:
\verb{||PE(tau)||_2^2 / (1-df(tau)/T)^2}
where
\verb{PE(tau) = E[R] - beta_\{S(tau)\} * OTFRP(tau)}
are the pricing errors of the model for given tuning parameter \code{tau},
with \code{S(tau)} being the index set of the nonzero Oracle TFRP computed with
tuning parameter \code{tau}, and
\verb{beta_\{S(tau)\} = Cov[R, F_\{S(tau)\}] * (Cov[F_\{S(tau)\}, R] * V[R]^\{-1\} * Cov[R, F_\{S(tau)\}])^\{-1\}}
the regression coefficients of the test assets excess returns on the
factor mimicking portfolios,
and \verb{df(tau) = |S(tau)|} are the degrees of freedom of the model, given by the
number of nonzero Oracle TFRP.
CV and RV, instead, choose the value of \code{tau} that minimize the criterium:
\verb{PE(tau)' * V[PE(tau)]^\{-1\} PE(tau)}
where \code{V[PE(tau)]} is the diagonal matrix collecting the marginal variances
of pricing errors \code{PE(tau)}, and each of these components are
aggregated over k-fold cross-validated data or over rolling windows of data,
respectively.
Oracle weights can be based on the correlation between factors and returns
(suggested approach),
on the regression coefficients of returns on factors or on the first-step
tradable risk premia estimator. Optionally computes the corresponding
heteroskedasticity and autocorrelation robust standard errors using the
Newey-West (1994) \url{doi:10.2307/2297912} plug-in procedure to select the number
of relevant lags, i.e., \code{n_lags = 4 * (n_observations/100)^(2/9)}.
For the standard error computations, the function allows to internally
pre-whiten the series by fitting a VAR(1),
i.e., a vector autoregressive model of order 1.
All details are found in Quaini-Trojani-Yuan (2023) \url{doi:10.2139/ssrn.4574683}.
}
\examples{
# import package data on 6 risk factors and 42 test asset excess returns
factors = intrinsicFRP::factors[,-1]
returns = intrinsicFRP::returns[,-1]

penalty_parameters = seq(0., 1., length.out = 100)

# compute optimal adaptive tradable factor risk premia and their standard errors
oracle_tfrp = OracleTFRP(
returns,
factors,
penalty_parameters,
include_standard_errors = TRUE
)

}
