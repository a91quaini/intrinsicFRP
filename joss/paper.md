---
title: 'intrinsicFRP: an R package for Oracle estimation and inference of tradable factor risk premia'
tags:
  - R
  - Factor Risk Premia
  - Testing Asset Pricing Models
  - Factor Selection
  - Oracle Estimation
authors:
  - name: Alberto Quaini
    orcid: 0000-0002-1251-0599
    affiliation: 1
    corresponding: true
affiliations:
 - name: Erasmus School of Economics, Erasmus University of Rotterdam
   index: 1
date: 18 September 2023
bibliography: paper.bib
---

# Summary
`intrinsicFRP` is an R package providing an efficient implementation of the
methods developed in [@quaini2023tradable] for the estimation and testing
of asset pricing models based on the notion of tradable factor risk premia.
Tradable factor risk premia are defined as the negative factor covariance with the
Stochastic Discount Factor projected on returns. Theoretically, they possess many 
desirable properties, namely they are well-identified even in the presence
of factors that are weakly correlated with returns (so-called weak factors), 
robust to potential misspecification of the underlying asset pricing model, 
and invariant to simple repackaging of the test 
assets. A simple estimator of tradable factor risk premia, obtained by replacing 
population moments with their empirical counterparts, is consistent and 
asymptotically Gaussian, even for weak factors having a correlation with returns 
vanishing at the rate $1/T^\alpha$ with $\alpha>1/2$, where $T$ is the sample size. 
For slowly-vanishing weak factors, i.e., with rate $1/T^{1/2}$, it's asymptotic distribution 
displays a bias term, which does not propagate to the other factor risk premia,
and it does not impair the overall asymptotic mean squared error of the estimator.
[@quaini2023tradable] further develop an Oracle estimator, which improves on the
simple estimator by employing an adaptively weighted lasso penalization to 
consistently remove the problematic weak factors, and it performs as well as if 
the weak or useless factors were known. More precisely, it consistently removes
the risk premium of any weak factor,
and it possess an efficient asymptotic Gaussian distribution for the remaining factors.

The aforementioned properties possessed by the simple and the Oracle tradable factor risk premia
estimators are remarkable, especially when compared to those of the commonly used
risk premia estimators in the literature, which are based on
usual notions of risk premia based on the Stochastic Discount Factor
projected on the factors; see, e.g., the Fama-MacBeth [@fama1973risk] and 
the misspecification-robust [@kan2013pricing] factor risk premia estimators.
In presence of weak factors, these estimators may be characterized by non-Gaussian
asymptotic distributions, and they may not even converge to zero for the weak factors.
Gaussian asymptotic distributions, and in some settings may even not be consistent; see, e.g., 
[@kan1999two], [@kleibergen2009tests] and [@gospodinov2014misspecification]. 

`intrinsicFRP` implements the methods for computing point estimates and standard 
errors for the simple and the Oracle tradable factor risk premia estimators,
as well as those of the commonly used Fama-MacBeth [@fama1973risk] and 
misspecification-robust [@kan2013pricing] factor risk premia. In order to 
maximize computing performances, all the methods are implemented in c++ and
rely on the [Armadillo](https://arma.sourceforge.net/) library for fast linear algebra calculations.
However, all methods are accessible from R functions, for ease of utilization.







`intrinsicFRP` is an R package Tradable factor risk premia are given by the negative factor covariance 
with the Stochastic Discount Factor projection on returns. This package 
provides efficient computation of tradable and Oracle tradable factor risk 
premia estimators and their standard errors; see A. Quaini, F. Trojani and 
M. Yuan (2023) <https://papers.ssrn.com/sol3/papers.cfm?abstract_id=4574683>. 
Tradable factor risk premia are robust to misspecification or weak 
identification in asset pricing models, and they are zero for any factor 
weakly correlated with returns. Their Oracle estimator performs as well as 
if the weak or useless factors were known in advance. This means it not only 
consistently removes useless factors and factors weakly correlated with 
returns but also gives rise to reliable tests of asset pricing models.
  [@quaini2023tradable] [@kan2013pricing] [fama1973risk]

# Statement of need

# Features

The R package `intrinsicFRP` comes with the following two datasets, which 
are sourced from the [Kenneth French data library](https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html):

- `returns`: monthly observations from `01/1970` to `12/2021` of the excess returns on the 25 Size/Book-to-Market and the 17 industry portfolios.
- `factors`:  monthly observations from `01/1970` to `12/2021` of the Fama-French 5 factors and the momentum factor.

The main functions in `intrinsicFRP` are:

- `TFRP()`: Computes tradable factor risk premia from data on factors and test asset excess returns. Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors.
- `OracleTFRP()`: Computes optimal Oracle tradable factor risk premia from data on factors and test asset excess returns. Tuning is performed via Generalized Cross Validation (GCV), Cross Validation (CV) or Rolling Validation (RV). Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors.
- `FRP()`: Computes Fama MachBeth (1973) or misspecification-robust factor risk premia of Kan Robotti Shanken (2013) from data on factors and test asset excess returns. Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors.
- `ChenFang2019BetaRankTest()`: Computes an asset pricing model identification test from data on factors and test asset excess returns. It employs a combination
of the rank tests in [@chen2019improved] and [@kleibergen2006generalized] on the matrix of regression coefficients of test asset excess returns on risk factors.
- `HJMisspecificationTest()`: Computes an asset pricing model misspecification test from data on factors and test asset excess returns. It employs the Hansen-Jagannatan test; see [@hansen1997assessing].

# Acknowledgments
The author would like to thank Fabio Trojani and Ming Yuan, co-creators 
together with the author of the methodologies underlying the R library `intrinsicFRP`.


# References
