---
title: 'intrinsicFRP: an R Package for Efficient Estimation and Inference of Asset Pricing Models'
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

intrinsicFRP is a lightweight R package providing an efficient implementation of methods from [@quaini2023tradable] for estimating and testing asset pricing models based on tradable factor risk premia. These premia are given by the negative factor covariance with the Stochastic Discount Factor projected onto returns. They are robust to model misspecification, asset repackaging and to the presence of factors that are weakly correlated with returns (weak factors).

A simple estimator for tradable factor risk premia, replacing population moments with empirical ones, remains consistent and asymptotically Gaussian, even for weak factors with correlations vanishing at a rate of $1/T^\alpha$, where $T$ is the sample size (for $\alpha>1/2$). For slowly vanishing weak factors ($1/T^{1/2}$), it exhibits a bias term, but this doesn't affect other factor risk premia.

[@quaini2023tradable] further introduces an Oracle estimator that enhances the simple one by using adaptively weighted lasso penalization to consistently remove problematic weak factors. It consistently removes the risk premium of any weak factor and possesses an efficient asymptotic Gaussian distribution for the remaining factors.

These properties of the simple and Oracle tradable factor risk premia estimators are noteworthy compared to commonly used estimators in the literature, like the Fama-MacBeth [@fama1973risk] and misspecification-robust [@kan2013pricing] estimators. These latter estimators, even in standard settings, may exhibit non-Gaussian asymptotic distributions and may not converge to zero for weak factors, as shown by in [@kan1999two], [@kleibergen2009tests], and [@gospodinov2014misspecification].

An illustrative figure below shows point estimates and standard errors for misspecification-robust (KRS), simple tradable (TFRP), and Oracle (O-TFRP) factor risk premia estimators. The computations use monthly observations from 01/1970 to 12/2021 for excess returns on 25 Size/Book-to-Market and 17 industry portfolios as test assets, along with the Fama-French 3 factors. An artificial useless factor is included, simulated from a standard normal distribution independent of asset returns.

In the figure, the misspecification-robust risk premium on the useless factor has the highest absolute value, with a wide 95% confidence interval. This reflects the poor statistical properties of the misspecification-robust estimator in the presence of weak factors. Conversely, the point estimate of the tradable risk premium on the useless factor remains close to zero, with a more concentrated 95% confidence interval. Finally, the point estimate of the Oracle tradable risk premium estimator on the useless factor is exactly zero, showcasing its variable selection capabilities.

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

The intrinsicFRP R package includes two datasets sourced from the Kenneth French data library [@french2012kenneth]:

- `returns`: Monthly observations from January 1970 to December 2021, containing excess returns data for 25 Size/Book-to-Market portfolios and 17 industry portfolios.
- `factors`: Monthly observations from January 1970 to December 2021, containing data for the Fama-French 5 factors and the momentum factor.

Additionally, intrinsicFRP implements methods for computing point estimates and standard errors for the simple and Oracle tradable factor risk premia estimators, as well as the Fama-MacBeth and misspecification-robust factor risk premia. The package also provides functions for conducting model identification tests developed in [@kleibergen2006generalized] and their bootstrap improvement in [@chen2019improved], along with the model misspecification test from [@hansen1997assessing]. The key functions in intrinsicFRP include:

- `TFRP()`: Computes tradable factor risk premia from data on factors and test asset excess returns. Optionally Optionally it computes robust standard errors accounting for heteroskedasticity and autocorrelation.
- `OracleTFRP()`: Computes optimal Oracle tradable factor risk premia using data on factors and test asset excess returns. Tuning is performed through Generalized Cross Validation (GCV), Cross Validation (CV), or Rolling Validation (RV). Robust standard errors can also be computed.
- `FRP()`: Computes Fama MachBeth (1973) or misspecification-robust factor risk premia following Kan Robotti Shanken (2013) based on data from factors and test asset excess returns. Robust standard errors can also be computed.
- `ChenFang2019BetaRankTest()`: Conducts an asset pricing model identification test using data on factors and test asset excess returns. This test combines rank tests from [@chen2019improved] and [@kleibergen2006generalized] on the matrix of regression coefficients of test asset excess returns against risk factors.
- `HJMisspecificationTest()`: Performs an asset pricing model misspecification test using data on factors and test asset excess returns. It utilizes the Hansen-Jagannatan test as described in [@hansen1997assessing].

To optimize computational performance, all these methods are implemented in C++ and make use of the Armadillo [@sanderson2016armadillo] library for efficient linear algebra calculations. However, for user convenience, the package's interface is entirely implemented in R.

# Dependencies

The `intrinsicFRP` library is a lightweight package with minimal dependencies, including:

- `Rcpp` [@eddelbuettel2018extending] and `RcppArmadillo` [@eddelbuettel2014rcpparmadillo]: They facilitate seamless integration between R, C++, and the armadillo C++ library.
- `graphics`: Provides R functions for creating basic graphics.
- `stats`: Offers R functions for performing statistical calculations and random number generation.

# Acknowledgments
The author extends gratitude to Fabio Trojani and Ming Yuan, who, along with the author, developed the methodologies underpinning the `intrinsicFRP` R library.


# References
