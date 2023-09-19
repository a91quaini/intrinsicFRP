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
`intrinsicFRP` is a lightweight R package providing an efficient implementation of the
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

A simple illustration of the properties discussed above is provided by the figure below, which reports the point estimates and standard errors 
of the misspecification-robust (KRS), tradable (TFRP) and Oracle (O-TFRP) factor risk premia estimators. computations are carried by considering
monthly observations from `01/1970` to `12/2021` of the
excess returns on the 25 Size/Book-to-Market and the 17 industry portfolios as the set of test assets, and the Fama-French 3 factors. The factor space is then augmented by means of an artificially constructed useless factor, which is simulated from a standard normal distribution and is independent of asset returns. 
These data are available in the package and are sourced from the Kenneth French data library [@french2012kenneth].

As shown in the figure, the misspecification-robust risk premium on the
useless factor has the highest absolute value, and its 95% confidence interval is very large. This reflects the poor statistical properties
that the misspecification-robust estimator has in the presence of weak factors. Instead, the point estimate of the tradable risk premium on the useless factor remains close to zero and its 95% confidence interval is much more concentrated than that of the misspecification-robust estimator.
Finally, the point estimate of the Oracle tradable risk premium estimator
on the useless factor is exactly zero. This reflects the Oracle variable selection properties possessed by the estimator. 

<p float="left">
<img src="../inst/examples/risk_premia.png" width="600" />
</p>
***Figure 1: **Estimates and 95% confidence intervals of the misspecification-robust (KRS), tradable (TFRP) and Oracle (O-TFRP) factor risk premia estimators for the market excess return (Mkt-RF), SMB, HML Fama-French factors and the simulated useless factor.*

# Statement of need

The R package `intrinsicFRP` is helpful for researchers and practitioners in finance due to the following pressing issues:

- Lack of Existing Code for Cutting-Edge Methods: The field of financial econometrics rapidly advancing. However, one major challenge is the absence of readily available code, especially in R, to implement the state-of-the-art methods. The `intrinsicFRP` package directly addresses this gap by providing a comprehensive and efficient implementation of the novel techniques established in [@quaini2023tradable].

- Dependency on Closed-Source Software: Historically, legacy code for factor models in finance has predominantly been written in proprietary programming languages like MATLAB. This dependence on closed-source software poses significant limitations, hindering accessibility, transparency, and collaboration within the research community. The `intrinsicFRP` R package is an open-source alternative that allows researchers to perform factor risk premia estimation and testing in a familiar, open environment.

By addressing these two key issues, the `intrinsicFRP` R package empowers researchers and practitioners in finance to:

- Readily apply the methods developed in [@quaini2023tradable] to their own research, facilitating the exploration of new avenues in asset pricing models and risk premia estimation. 

- Promote open-source science software and overcome MATLAB dependencies, thereby improving transparency and reproducibility of research in the Asset Pricing and related fields. By eliminating software constraints, we aim to facilitate access to advanced methods for finance researchers and practitioners, fostering collaboration among researchers and enabling the wider financial community to further develop and validate the suggested approaches.

# Features

The R package `intrinsicFRP` comes with the following two datasets, which 
are sourced from the Kenneth French data library [@french2012kenneth]:

- `returns`: monthly observations from `01/1970` to `12/2021` of the excess returns on the 25 Size/Book-to-Market and the 17 industry portfolios.
- `factors`:  monthly observations from `01/1970` to `12/2021` of the Fama-French 5 factors and the momentum factor.

Moreover, `intrinsicFRP` implements the methods for computing point estimates and standard 
errors for the simple and the Oracle tradable factor risk premia estimators,
as well as those of the Fama-MacBeth and 
misspecification-robust factor risk premia. The package further 
implements the rank test for model identification developed in 
[@kleibergen2006generalized] and its bootstrap improvement in [@chen2019improved],
as well as the model misspecification test of [@hansen1997assessing]. The main functions in `intrinsicFRP` are:

- `TFRP()`: Computes tradable factor risk premia from data on factors and test asset excess returns. Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors.
- `OracleTFRP()`: Computes optimal Oracle tradable factor risk premia from data on factors and test asset excess returns. Tuning is performed via Generalized Cross Validation (GCV), Cross Validation (CV) or Rolling Validation (RV). Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors.
- `FRP()`: Computes Fama MachBeth (1973) or misspecification-robust factor risk premia of Kan Robotti Shanken (2013) from data on factors and test asset excess returns. Optionally computes the corresponding heteroskedasticity and autocorrelation robust standard errors.
- `ChenFang2019BetaRankTest()`: Computes an asset pricing model identification test from data on factors and test asset excess returns. It employs a combination
of the rank tests in [@chen2019improved] and [@kleibergen2006generalized] on the matrix of regression coefficients of test asset excess returns on risk factors.
- `HJMisspecificationTest()`: Computes an asset pricing model misspecification test from data on factors and test asset excess returns. It employs the Hansen-Jagannatan test; see [@hansen1997assessing].

In order to 
maximize computing performances, all the methods are implemented in C++ and
rely on the Armadillo [@sanderson2016armadillo] library for fast linear algebra calculations.
However, for ease of utilization, the user interface is completely implemented in R.

# Dependencies

`intrinsicFRP` is a lightweight library only depending on:

- `Rcpp` [@eddelbuettel2018extending] and `RcppArmadillo` [@eddelbuettel2014rcpparmadillo]: 
packages for seamless integration of R, C++ and the `armadillo` C++ library.
- `graphics`: R functions for base graphics.
- `stats` : R functions for statistical calculations and random number generation.

# Acknowledgments
The author would like to thank Fabio Trojani and Ming Yuan, co-creators 
together with the author of the methodologies underlying the R library `intrinsicFRP`.


# References
