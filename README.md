---
bibliography: readme.bib
---

# intrinsicFRP: An R Package for Factor Model Asset Pricing

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/intrinsicFRP)](https://CRAN.R-project.org/package=intrinsicFRP)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/a91quaini/intrinsicFRP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/a91quaini/intrinsicFRP/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/a91quaini/intrinsicFRP/branch/main/graph/badge.svg?token=0F8R40B0FP)](https://app.codecov.io/gh/a91quaini/intrinsicFRP)
[![](https://cranlogs.r-pkg.org/badges/grand-total/intrinsicFRP)](https://cran.rstudio.com/web/packages/intrinsicFRP/index.html)
<!-- badges: end -->

Author: Alberto Quaini

The `intrinsicFRP` library implements **functions designed for comprehensive evaluation and testing of asset pricing models**, focusing on the estimation and assessment of factor risk premia, selection of "useful" risk factors (those displaying non-zero population correlation with test asset returns), heteroskedasticity and autocorrelation robust covariance matrix estimation, examination of model misspecification, and validation of model identification.

For **estimating and testing factor risk premia**, our toolkit incorporates the classic two-pass approach outlined by [@fama1973risk], the misspecification-robust methodology proposed by [@kan2013pricing], and the tradable and "Oracle" tradable approaches introduced by [@quaini2023tradable].

The "Oracle" estimator of [@quaini2023tradable] can also be used to **identify the set of useful risk factors**. The motive is that, 
using the terminology of [@fan2001variable],
this estimator has the so-called "Oracle" variable selection property, i.e., 
it consistently selects the useful factors. More precisely, the
probability that the factors selected by the estimator are indeed useful factors
tends to 1 as the sample size tends to infinity.

For **heteroskedasticity and autocorrelation robust covariance matrix estimation**, 
the package implements the [@newey1994automatic] covariance estimator.

For **evaluating model misspecification**, the toolkit implements the HJ model misspecification test formulated by [@hansen1997assessing], following the implementation of [@kan2008specification].

Lastly, the functions for **testing model identification** are specialized versions of the rank tests proposed by [@kleibergen2006generalized] and [@chen2019improved]. These tests are specifically tailored to assess the regression coefficient matrix of test asset returns on risk factors.
  
## Citation

To cite `intrinsicFRP` in publications, please use:

> Quaini, A. (2023). `intrinsicFRP`: R Package For Factor Model Asset Pricing. `R` package version 1.0.0. URL: https://cran.r-project.org/web/packages/intrinsicFRP/index.html.

## Installation

### Building from CRAN

Package `intrinsicFRP` is on CRAN (The Comprehensive R Archive Network),
hence the latest release can be easily installed from the `R` command
line via

```R
install.packages("intrinsicFRP")
```

### Building from source

To install the latest (possibly unstable) development version from
GitHub, you can pull this repository and install it from the `R` command
line via

```R
# if you already have package `devtools` installed, you can skip the next line
install.packages("devtools")
devtools::install_github("a91quaini/intrinsicFRP")
```

Package `intrinsicFRP` contains `C++` code that needs to be
compiled, so you may need to download and install the [necessary tools
for MacOS](https://cran.r-project.org/bin/macosx/tools/) or the
[necessary tools for
Windows](https://cran.r-project.org/bin/windows/Rtools/).

## Features

### Functions

R package `intrinsicFRP` implements the following functions:

- `FRP()`: Computes the Fama-MachBeth (1973) [<doi:10.1086/260061>] factor risk premia: `FMFRP = (beta' * beta)^{-1} * beta' * E[R]` where `beta = Cov[R, F] * V[F]^{-1}` or the misspecification-robust factor risk premia of [@kan2008specification]: `KRSFRP = (beta' * V[R]^{-1} * beta)^{-1} * beta' * V[R]^{-1} * E[R]` from data on factors `F` and test asset excess returns `R`. Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors using the [@newey1994automatic] plug-in procedure to select the number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
- `TFRP()`: Computes tradable factor risk premia from data on factors `F` and test asset excess returns `R`: `TFRP = Cov[F, R] * Var[R]^{-1} * E[R]`. Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors using the [@newey1994automatic] plug-in procedure to select the number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`. All details are found in [@quaini2023tradable].
- `OracleTFRP()`: Computes optimal adaptive tradable factor risk premia for various penalty parameter values `tau` from data on `K` factors `F = [F_1,...,F_K]'` and test asset excess returns `R`: `OTFRP = argmin_x ||TFRP - x||_2^2 + tau * sum_{k=1}^K w_k * |x_k|` where `TFRP` is the sample tradable factor risk premia estimator  and the "Oracle" weights are given by `w_k = 1/||corr[F_k, R]||_2^2`. Tuning of the penalty parameter is performed via Generalized Cross Validation (GCV), Cross Validation (CV) or Rolling Validation (RV).
Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors using the [@newey1994automatic] plug-in procedure to select the number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`. All details are found in Quaini-Trojani-Yuan (2023) [<doi:10.2139/ssrn.4574683>].
- `HACcovariance()`: estimates the long-run covariance matrix of a multivariate centred time series accounting for heteroskedasticity and autocorrelation using the Newey-West estimator. If the number of lags is not provided, they are selected using the [@newey1994automatic] plug-in procedure, where `n_lags = 4 * (n_observations/100)^(2/9)`. The function allows to internally prewhiten the series by fitting a VAR(1).
- `HJMisspecificationTest()`: Computes the [@hansen1997assessing] model misspecification statistic and p-value of an asset pricing model from test asset excess returns `R` and risk factors `F`: `HJDISTANCE = min_{d} (E[R] - beta * d)' * Var[R]^{-1} * (E[R] - beta * d)` where `beta = Cov[R, F] * Var[F]^{-1}` are the regression coefficients of test asset excess returns `R` on risk factors `F`. Details on the closed-form solution and computation of the p-values are found in [@kan2008specification].
- `ChenFang2019BetaRankTest()`: Computes the test statistic and p-value of the [@chen2019improved] rank test of the null that the matrix of regression loadings of `N` test asset excess returns `R` on `K` risk factors `F`, with `N > K`, has reduced rank: `H: rank(beta) < K` where `beta = Cov[R, F] * Var[F]^{-1}` is the `N x K` matrix of regression coefficients. If `target_level_kp2006_rank_test > 0`, it uses the iterative [@kleibergen2009tests] rank test to estimate the initial rank, with `level = target_level_kp2006_rank_test / n_factors`. If `target_level_kp2006_rank_test <= 0`, the initial rank estimator is taken to be the number of singular values above `n_observations^(-1/4)`.
- `IterativeKleibergenPaap2006BetaRankTest()`: Computes the test statistic and p-values of the iterative [@kleibergen2009tests] rank test of the null that the matrix of regression loadings of `N` test asset excess returns `R` on `K` risk factors `F`, with `N > K`, has reduced rank: `H: rank(beta) < K` where `beta = Cov[R, F] * Var[F]^{-1}` is the `N x K` matrix of regression coefficients. Since it is an iterative procedure, testing `H: rank(beta) = q` with `q=0,...,K-1`, it internally employs a Bonferroni correction for multiple testing. It also returns an estimate of the rank, computed as the first value `q` with associated p-value below the given `level = target_level / n_factors`.
  
For usage details, type `?FunctionName` in the R console, e.g.:

```R
?TFRP
```

## Data

The `intrinsicFRP` R package includes a dataset comprising following
test asset excess returns and risk factors frequently used 
in the asset pricing literature:

- `returns`: Monthly observations from January 1970 to December 2021, containing excess returns data for 25 Size/Book-to-Market portfolios and 17 industry portfolios.
- `factors`: Monthly observations from January 1970 to December 2021, containing data for the Fama-French 5 factors and the momentum factor.

This dataset was sourced from the [Kenneth French data library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html)
and processed so that the observations on the test assets and the factors are not
expressed in percentage points, and that returns on the test assets are in excess 
of the risk-free rate.

## Examples

## Real data example 1: estimation and inference of factor risk premia

Compute various factor risk premia estimates and corresponding standard errors for the Fama-French 6 factors and a (simulated) "useless" factor.

```R
# import package data on 6 risk factors and 42 test asset excess returns
# remove the first column containing the date
factors = intrinsicFRP::factors[,-1]
returns = intrinsicFRP::returns[,-1]

# simulate a useless factor and add it to the matrix of factors
set.seed(23)
factors = cbind(
  factors,
  stats::rnorm(n = nrow(factors), sd = stats::sd(factors[,3]))
)
colnames(factors) = c(colnames(intrinsicFRP::factors[,1:6]), "Useless")

# index set of specific factor models
# Fama-French 3 factor model
ff3 = 1:3
# Fama-French 6 factor model
ff6 = 1:6         # "Mkt-RF" "SMB" "HML" "RMW" "CMA" "Mom"
# model comprising the Fama-French 6 factors and the simulated useless factor
ff6usl = 1:7      # "Mkt-RF" "SMB" "HML" "RMW" "CMA" "Mom" "Useless"

# compute tradable factor risk premia and their standard errors
tfrp = intrinsicFRP::TFRP(returns, factors[,ff6usl], include_standard_errors = TRUE)

# compute the GLS factor risk premia of Kan Robotti and Shanken (2013) and their
# standard errors
krs_frp = intrinsicFRP::FRP(returns, factors[,ff6usl], include_standard_errors = TRUE)

# set penalty parameters
penalty_parameters = seq(1e-4, 1e-2, length.out = 1000)

# compute Oracle tradable factor risk premia and their standard errors
# for low factor models, no need for the "one standard deviation" tuning rule
oracle_tfrp = intrinsicFRP::OracleTFRP(
  returns,
  factors[,ff6usl],
  penalty_parameters,
  include_standard_errors = TRUE,
  one_stddev_rule = FALSE
)
```

<!--```R
# create dataframe
df <- data.frame(
  Factor = factor(
    rep(colnames(factors[,ff6usl]), 3),
    levels = colnames(factors[,ff6usl])
  ),
  Estimator = factor(
    rep(c("KRS-FRP", "TFRP", "O-TFRP"), each=ncol(factors[,ff6usl])),
    levels = c("KRS-FRP", "TFRP", "O-TFRP")
  ),
  risk_premia = c(krs_frp$risk_premia, tfrp$risk_premia, oracle_tfrp$risk_premia),
  standard_errors = c(
    krs_frp$standard_errors, tfrp$standard_errors, oracle_tfrp$standard_errors
  )
)

# Create the plot
ggplot2::ggplot(df, ggplot2::aes(
  x = as.factor(.data$Factor), y = .data$risk_premia, fill = .data$Estimator)) +
  ggplot2::theme(text=ggplot2::element_text(size=16)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge", width=0.5, color="black") +
  ggplot2::labs(x = "Factor", y = "Risk Premia") +
  ggplot2::geom_errorbar(ggplot2::aes(
    x=as.factor(Factor),
    ymin=risk_premia - stats::qnorm(0.975) * standard_errors,
    ymax=risk_premia + stats::qnorm(0.975) * standard_errors),
    linewidth=.8, position = ggplot2::position_dodge(0.5), width = 0.25)
```-->

Tuning model score of the Oracle TFRP estimator:

<p float="left">
<img src="inst/examples/model_score.png" width="600" />
</p>

Visualization of the misspecification-robust factor risk premia (KRS-FRP), tradable factor risk premia (TFRP) and Oracle TFRP (O-TFRP) estimates:

<p float="left">
<img src="inst/examples/risk_premia.png" width="600" />
</p>

## Real data example 2: testing misspecification and identification of asset pricing models

Compute the HJ misspecification test of the Fama-French 6 factor model and identification tests of the Fama-French 6 factor model and the (unidentified)
model comprising the Fama-French 6 factors and the simulated useless factor.

```R
# compute the HJ misspecification test of the Fama-French 3 and 6 factor models
intrinsicFRP::HJMisspecificationTest(returns, factors[,ff3])["p-value"]
intrinsicFRP::HJMisspecificationTest(returns, factors[,ff6])["p-value"]

# compute identification tests of the Fama-French 6 factor model
intrinsicFRP::IterativeKleibergenPaap2006BetaRankTest(returns, factors[,ff6])
intrinsicFRP::ChenFang2019BetaRankTest(returns, factors[,ff6])["p-value"]

# compute identification tests of unidentified factor model comprising the
# Fama-French 6 factors and the simulated useless factor
intrinsicFRP::IterativeKleibergenPaap2006BetaRankTest(returns, factors[,ff6usl])
intrinsicFRP::ChenFang2019BetaRankTest(returns, factors[,ff6usl])["p-value"]
```

Result of the HJ misspecification test:
```R
# HJ misspecification test p-value for the Fama-French 3 factor model
1.526582e-07

# HJ misspecification test p-value for the Fama-French 6 factor model
2.133855e-05
```

Since the p-value of both HJ misspecification tests is below the standard thresholds of $10\%$, $5\%$ and $1\%$, we reject the Null that the Fama-French 3 and 6 factor models are correctly specified.

Results of the identification tests for the Fama-French 6 factor model:
```R
# output of the Iteraive Kleibergen Paap (2006) Beta Rank Test
$rank
[1] 6

$q
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0    1    2    3    4    5

$statistics
[1] 425512.0961 156257.6897  36635.6923   1955.2522    486.3576    191.7194

$pvalues
[1] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 3.986757e-61 9.076916e-23

# p-value of the Chen Fang (2019) Beta Rank Test
$`p-value`
[1] 0
```

Since the largest p-value of the Iteraive Kleibergen Paap (2006) Beta Rank Test and the p-value of the Chen Fang (2019) Beta Rank Test are below the standard thresholds of $10\%$, $5\%$ and $1\%$, we reject the Null that the Fama-French 6 factor model is not identified.

For sanity check, compute identification tests of the unidentified model comprising the Fama-French 6 factors and a (simulated) "useless" factor.

```R
# compute identification test of unidentified factor models
intrinsicFRP::IterativeKleibergenPaap2006BetaRankTest(returns, factors[,ff6usl])
intrinsicFRP::ChenFang2019BetaRankTest(returns, factors[,ff6usl])["p-value"]
```

Results:
```R
# output of the Iteraive Kleibergen Paap (2006) Beta Rank Test
$rank
[1] 6

$q
     [,1] [,2] [,3] [,4] [,5] [,6] [,7]
[1,]    0    1    2    3    4    5    6

$statistics
[1] 505705.19553 171847.08238  42170.13669   2214.27454    591.29522    274.41845
[7]     58.45019

$pvalues
[1] 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 1.600124e-65 8.255293e-25
[7] 1.039242e-02

# p-value of the Chen Fang (2019) Beta Rank Test
$`p-value`
[1] 0.166
```

Since the largest p-value of the Iteraive Kleibergen Paap (2006) Beta Rank Test and the p-value of the Chen Fang (2019) Beta Rank Test are above the standard thresholds of $10\%$, $5\%$ and $1\%$, we do not reject the Null that the Fama-French 6 factor model is not identified.

## Dependencies

To optimize computational performance, all methods implemented in package `intrinsicFRP` are written in C++ and make use of the [Armadillo](https://arma.sourceforge.net/) library for efficient linear algebra calculations. However, for user convenience, the interface of package `intrinsicFRP` is entirely implemented in R, with minimal dependencies, including:

- `Rcpp` [@eddelbuettel2018extending] and `RcppArmadillo` [@eddelbuettel2014rcpparmadillo]: They facilitate seamless integration between R, C++, and the armadillo C++ library.
- `graphics`: Provides R functions for creating basic graphics.
- `stats`: Offers R functions for performing statistical calculations and random number generation.

## Issues, bug reports, contributions, further help

You can raise issues, report bugs, seek for further help, or submit your contribution to the R package `intrinsicFRP` at the
github repository [github.com/a91quaini/intrinsicFRP](github.com/a91quaini/intrinsicFRP). 

For bug reports, you are kindly asked to make a small and self-contained program which exposes the bug.

## References
