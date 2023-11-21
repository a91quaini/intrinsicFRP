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

ggplot2::ggsave(
  "inst/examples/risk_premia.png",
  width = 7,
  height = 5,
  dpi=600
)

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
