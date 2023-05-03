# intrinsicFRP

<!-- badges: start -->
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!--[![codecov](https://codecov.io/gh/a91quaini/adaptiveIFRP/branch/main/graph/badge.svg?token=0F8R40B0FP)](https://codecov.io/gh/a91quaini/adaptiveIFRP)-->
[![R-CMD-check](https://github.com/a91quaini/intrinsicFRP/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/a91quaini/intrinsicFRP/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Author: Alberto Quaini

Efficient implementation of intrinsic factor risk premia and their
“adaptive” version based on the methods developed in Quaini, Trojani and
Yuan 2023 “Intrinisic Factor Risk Premia and Testing of Asset Pricing
Models”.

## Installation

You can install the development version of intrinsicFRP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("a91quaini/intrinsicFRP")
```

## Example

Compute risk premia estimates and corresponding standard errors:

```R
library(intrinsicFRP)

# import package data on 15 risk factors and 42 test asset excess returns
factors = factors[,-1]
returns = returns[,-1]

# compute intrinsic factor risk premia and their standard errors
ifrp = IFRP(returns, factors, include_standard_errors = TRUE)

# set penalty parameters
penalty_parameters = seq(0., .01, length.out = 100)

# compute optimal adaptive intrinsic factor risk premia and their standard
# errors
aifrp = OptimalAdaptiveIFRP(
  returns,
  factors,
  penalty_parameters,
  include_standard_errors = TRUE
)

# compute the GLS factor risk premia of Kan Robotti and Shanken (2013) and their
# standard errors
krs_frp = FRP(returns, factors, include_standard_errors = TRUE)
```

Output:

<!--<p float="left">
<img src="inst/examples/risk_premia.svg" width="800" />
</p>-->

## Functions

For usage details, type ?FunctionName in the R console, e.g.:

```R
?IFRP
```

### Main functions:

- **IFRP**: Computes intrinsic factor risk premia from data on factors
  and test asset excess returns. Optionally computes the 
  correspodning heteroskedasticity and autocorrelation robust standard errors 
  using the Newey-West (1994) plug-in procedure to select the number of 
  relevant lags, i.e., n_lags = 4 * (n_observations/100)^(2/9).
- **OptimalAdaptiveIFRP**: Computes optimal adaptive intrinsic factor risk 
  premia based on pre-computed intrinsic factor risk premia and adaptive penalty 
  weights for various penalty parameter values. Tuning is performed via
  Generalized Cross Validation (GCV), Cross Validation (CV) or Rolling 
  Validation (RV). Adaptive weights are based on the correlation between factors 
  and returns, on the regression coefficients of returns on factors or on the 
  first-step intrinsic risk premia estimator. Optionally computes the 
  corresponding heteroskedasticity and autocorrelation robust standard errors 
  using the Newey-West (1994) plug-in procedure to select the number of 
  relevant lags, i.e., n_lags = 4 * (n_observations/100)^(2/9).
- **AdaptiveIFRP**: Computes adaptive intrinsic factor risk premia based on
  pre-computed intrinsic factor risk premia and adaptive penalty weights for
  various penalty parameter values.
- **FRP**: Computes Fama-MachBeth or misspecification-robust factor risk
  premia of Kan Robotti Shanken (2013) from data on factors and test
  asset excess returns. Optionally computes the 
  corresponding heteroskedasticity and autocorrelation robust standard errors 
  using the Newey-West (1994) plug-in procedure to select the number of 
  relevant lags, i.e., n_lags = 4 * (n_observations/100)^(2/9).
- **PlotAdaptiveIFRPModelScore**: Plots the model score of the adaptive IFRP for
  each penalty parameter value, highlighting the minimum attained model score 
  and the optimal one (which differs from the minimum in case the optimal 
  parameter was computed with the `one_stddev_rule` option TRUE).

## References

Quaini, A., Trojani, F. and Yuan, M., 2023. Intrinsic Factor Risk Premia
and Testing of Asset Pricing Models. Working paper.

Kan, R., Robotti, C. and Shanken, J., 2013. Pricing model performance and the two‐pass cross‐sectional regression methodology. The Journal of Finance, 68(6), pp.2617-2649.
