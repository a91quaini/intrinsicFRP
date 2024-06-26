# if you already have package `stats` and `ggplot2` installed
# you can skip the next line
install.packages(c("stats", "ggplot2"))

# import package data on 6 risk factors and 42 test asset excess returns
# remove the first column containing the date
factors = intrinsicFRP::factors[,-1]
returns = intrinsicFRP::returns[,-1]
RF = intrinsicFRP::risk_free[,-1]

# simulate a useless factor and add it to the matrix of factors
set.seed(23)
factors = cbind(
  factors,
  stats::rnorm(n = nrow(factors), sd = stats::sd(factors[,3]))
)
colnames(factors) = c(colnames(intrinsicFRP::factors[,2:7]), "Useless")

# index set of specific factor models
# Fama-French 3 factor model
ff3 = 1:3
# Fama-French 6 factor model
ff6 = 1:6         # "Mkt-RF" "SMB" "HML" "RMW" "CMA" "Mom"
# model comprising the Fama-French 6 factors and the simulated useless factor
ff6usl = 1:7      # "Mkt-RF" "SMB" "HML" "RMW" "CMA" "Mom" "Useless"

# compute the factor SDF coefficients and their standard errors
# for the Fama-MacBeth two-pass procedure
fm_sdf = intrinsicFRP::SDFCoefficients(
  returns,
  factors[,ff6usl],
  misspecification_robust = FALSE,
  include_standard_errors = TRUE
)
# and for the misspecification-robust procedure of Gospodinov Kan and Robottu
gkr_sdf = intrinsicFRP::SDFCoefficients(
  returns, factors[,ff6usl],
  include_standard_errors = TRUE
)

# create dataframe
df = data.frame(
  Factor = factor(
    rep(colnames(factors[,ff6usl]), 2),
    levels = colnames(factors[,ff6usl])
  ),
  Estimator = factor(
    rep(c("FM", "GKR"), each=ncol(factors[,ff6usl])),
    levels = c("FM", "GKR")
  ),
  sdf_coefficients = c(fm_sdf$sdf_coefficients, gkr_sdf$sdf_coefficients),
  standard_errors = c(fm_sdf$standard_errors, gkr_sdf$standard_errors)
)

# Create the plot
ggplot2::ggplot(df, ggplot2::aes(
  x = as.factor(.data$Factor), y = .data$sdf_coefficients, fill = .data$Estimator)) +
  ggplot2::theme(text=ggplot2::element_text(size=16)) +
  ggplot2::geom_bar(stat = "identity", position = "dodge", width=0.5, color="black") +
  ggplot2::labs(x = "Factor", y = "SDF coefficient") +
  ggplot2::geom_errorbar(ggplot2::aes(
    x=as.factor(Factor),
    ymin=sdf_coefficients - stats::qnorm(0.975) * standard_errors,
    ymax=sdf_coefficients + stats::qnorm(0.975) * standard_errors),
    linewidth=.8, position = ggplot2::position_dodge(0.5), width = 0.25)

ggplot2::ggsave(
  "inst/examples/sdf_coefficients.png",
  width = 7,
  height = 5,
  dpi=600
)


# compute tradable factor risk premia and their standard errors
tfrp = intrinsicFRP::TFRP(returns, factors[,ff6usl], include_standard_errors = TRUE)

# compute the Fama-MacBeth factor risk premia and their standard errors
fm_frp = intrinsicFRP::FRP(returns, factors[,ff6usl], include_standard_errors = TRUE, misspecification_robust = FALSE)

# compute the GLS factor risk premia of Kan Robotti and Shanken (2013) and their
# standard errors
krs_frp = intrinsicFRP::FRP(returns, factors[,ff6usl], include_standard_errors = TRUE)

# set penalty parameters
penalty_parameters = seq(1e-4, 4e-3, length.out = 1000)

# compute Oracle tradable factor risk premia and their standard errors
# for low factor models, no need for the "one standard deviation" tuning rule
oracle_tfrp = intrinsicFRP::OracleTFRP(
  returns,
  factors[,ff6usl],
  penalty_parameters,
  include_standard_errors = TRUE,
  one_stddev_rule = FALSE
)

# create dataframe
df = data.frame(
  Factor = factor(
    rep(colnames(factors[,ff6usl]), 4),
    levels = colnames(factors[,ff6usl])
  ),
  Estimator = factor(
    rep(c("FM", "KRS", "TFRP", "O-TFRP"), each=ncol(factors[,ff6usl])),
    levels = c("FM", "KRS", "TFRP", "O-TFRP")
  ),
  risk_premia = c(fm_frp$risk_premia, krs_frp$risk_premia, tfrp$risk_premia, oracle_tfrp$risk_premia),
  standard_errors = c(
    fm_frp$standard_errors, krs_frp$standard_errors, tfrp$standard_errors, oracle_tfrp$standard_errors
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

ggplot2::ggsave(
  "inst/examples/risk_premia.png",
  width = 9,
  height = 5,
  dpi=600
)

# recover the indices of the factors selected by the Oracle TFRP estimator
which(oracle_tfrp$risk_premia != 0)

# compute the GKR factor screening procedure
intrinsicFRP::GKRFactorScreening(returns, factors[,ff6])

# compute the FGX factor screening procedure
gross_returns = returns + 1 + RF
fgx_output = intrinsicFRP::FGXFactorsTest(gross_returns, factors[,ff3], factors[,4:7])

alpha = 0.05
quantile = stats::qnorm(1 - alpha / 2)
lows = fgx_output$sdf_coefficients - fgx_output$standard_errors * quantile
highs = fgx_output$sdf_coefficients + fgx_output$standard_errors * quantile
# recover the indices of the factors selected by the FGX procedure
which(lows * highs > 0)

# compute the HJ misspecification distance of the Fama-French 3 and 6 factor models
intrinsicFRP::HJMisspecificationDistance(returns, factors[,ff3])
intrinsicFRP::HJMisspecificationDistance(returns, factors[,ff6])

# compute identification tests of the Fama-French 6 factor model
intrinsicFRP::IterativeKleibergenPaap2006BetaRankTest(returns, factors[,ff6])
intrinsicFRP::ChenFang2019BetaRankTest(returns, factors[,ff6])

# compute identification tests of unidentified factor model comprising the
# Fama-French 6 factors and the simulated useless factor
intrinsicFRP::IterativeKleibergenPaap2006BetaRankTest(returns, factors[,ff6usl])
intrinsicFRP::ChenFang2019BetaRankTest(returns, factors[,ff6usl])

# compute the HAC covariance matrix of the residuals from a regression
# of returns on factors
residuals = stats::lm.fit(factors, returns)$residuals
asy_cov_residuals = intrinsicFRP::HACcovariance(residuals)
