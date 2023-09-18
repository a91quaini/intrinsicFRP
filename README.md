
# intrinsicFRP: Oracle tradable factor risk premia

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/intrinsicFRP)](https://CRAN.R-project.org/package=intrinsicFRP)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![R-CMD-check](https://github.com/a91quaini/intrinsicFRP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/a91quaini/intrinsicFRP/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/a91quaini/intrinsicFRP/branch/main/graph/badge.svg?token=0F8R40B0FP)](https://app.codecov.io/gh/a91quaini/intrinsicFRP)
<!-- badges: end -->

Author: Alberto Quaini

Efficient implementation of the sample and Oracle tradable factor risk premia estimators based on the methods developed in A. Quaini, F. Trojani and
M. Yuan (2023) [Tradable Factor Risk Premia and Oracle Tests of Asset Pricing
Models](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4574683).

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

### Main functions:

- `TFRP()` [Prevously: `IFRP`]: Computes tradable factor risk premia from data on factors
  and test asset excess returns. Optionally computes the 
  corresponding heteroskedasticity and autocorrelation robust standard errors 
  using the Newey-West (1994) plug-in procedure to select the number of 
  relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
- `OracleTFRP()` [Prevously: `OptimalAdaptiveIFRP`]: Computes optimal Oracle 
tradable factor risk 
  premia for various penalty parameter values from data on factors and test 
  asset excess returns. Tuning is performed via Generalized Cross Validation (GCV),
  Cross Validation (CV) or Rolling Validation (RV). Adaptive weights are based 
  on the correlation between factors and returns, on the regression coefficients 
  of returns on factors or on the first-step intrinsic risk premia estimator.
  Optionally computes the corresponding heteroskedasticity and autocorrelation 
  robust standard errors using the Newey-West (1994) plug-in procedure to select 
  the number of relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
- `FRP()`: Computes Fama MachBeth (1973) or misspecification-robust factor risk
  premia of Kan Robotti Shanken (2013) from data on factors and test
  asset excess returns. Optionally computes the 
  corresponding heteroskedasticity and autocorrelation robust standard errors 
  using the Newey-West (1994) plug-in procedure to select the number of 
  relevant lags, i.e., `n_lags = 4 * (n_observations/100)^(2/9)`.
- `ChenFang2019BetaRankTest()`: Computes the Chen fang 2019 rank statistic and p-value of the
null that the matrix of regression loadings of test asset excess returns on
risk factors has reduced rank.
- `HJMisspecificationTest()`: Computes the Hansen-Jagannatan misspecification statistic and p-value of an asset pricing model from test asset excess returns and risk factors.
  
For usage details, type `?FunctionName` in the R console, e.g.:

```R
?TFRP
```

## Example

Compute risk premia estimates and corresponding standard errors:

```R
# import package data on 6 risk factors and 42 test asset excess returns
factors = intrinsicFRP::factors[,c(2,3,4)]
returns = intrinsicFRP::returns[,-1]

# simulate a useless factor and add it to the matrix of factors
set.seed(23)
factors = cbind(factors, stats::rnorm(nrow(factors), sd = stats::sd(factors[,1])))
colnames(factors) = c(colnames(intrinsicFRP::factors[,c(2,3,4)]), "Usless")

# compute tradable factor risk premia and their standard errors
tfrp = TFRP(returns, factors, include_standard_errors = TRUE)

# compute the GLS factor risk premia of Kan Robotti and Shanken (2013) and their
# standard errors
krs_frp = FRP(returns, factors, include_standard_errors = TRUE)

# set penalty parameters
penalty_parameters = seq(1e-4, 1e-2, length.out = 1000)

# compute Oracle tradable factor risk premia and their standard errors
oracle_tfrp = OracleTFRP(
  returns,
  factors,
  penalty_parameters,
  include_standard_errors = TRUE,
)

```

<!--```R
# create dataframe
df <- data.frame(
  Factor = factor(
    rep(colnames(factors[,1:4]), 3),
    levels = colnames(factors[,1:4])
  ),
  Estimator = rep(c("TFRP", "O-TFRP", "KRS"), each=ncol(factors[,1:4])),
  risk_premia = c(tfrp$risk_premia, oracle_tfrp$risk_premia, krs_frp$risk_premia),
  standard_errors = c(
    tfrp$standard_errors, oracle_tfrp$standard_errors, krs_frp$standard_errors
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
Oracle TFRP model score:

<p float="left">
<img src="inst/examples/model_score.png" width="600" />
</p>

Compare risk premia estimates of misspecification-robust factor risk premia (KRS), tradable factor risk premia (TFRP) and Oracle TFRP (O-TFRP):

<p float="left">
<img src="inst/examples/risk_premia.png" width="600" />
</p>

## References

Quaini, A., Trojani, F. and Yuan, M., 2023. [Tradable Factor Risk Premia
and Oracle Tests of Asset Pricing Models](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4574683).

Kan, R., Robotti, C. and Shanken, J., 2013. Pricing model performance and the two‐pass cross‐sectional regression methodology. The Journal of Finance, 68(6), pp.2617-2649.

Fama, E.F. and MacBeth, J.D., 1973. Risk, return, and equilibrium: Empirical tests. Journal of political economy, 81(3), pp.607-636.
