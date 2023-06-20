# import package data on 6 risk factors and 42 test asset excess returns
factors = intrinsicFRP::factors[,c(2,3,4)]
returns = intrinsicFRP::returns[,-1]

# simulate a useless factor and add it to the matrix of factors
factors = cbind(factors, stats::rnorm(nrow(factors), sd = stats::sd(factors[,1])))
colnames(factors) = c(colnames(intrinsicFRP::factors[,c(2,3,4)]), "Usless")

# compute intrinsic factor risk premia and their standard errors
ifrp = IFRP(returns, factors, include_standard_errors = TRUE)

# compute the GLS factor risk premia of Kan Robotti and Shanken (2013) and their
# standard errors
krs_frp = FRP(returns, factors, include_standard_errors = TRUE)

# set penalty parameters
penalty_parameters = seq(1e-4, 1e-1, length.out = 1000)

# compute optimal adaptive intrinsic factor risk premia and their standard
# errors
aifrp = OptimalAdaptiveIFRP(
  returns,
  factors,
  penalty_parameters,
  include_standard_errors = FALSE,
)

afrp = OptimalAdaptiveFRP(
  returns,
  factors,
  penalty_parameters,
  include_standard_errors = TRUE
)

# create dataframe
df <- data.frame(
  Factor = rep(colnames(factors[,1:4]), 3),
  Estimator = rep(c("IFRP", "A-IFRP", "KRS"), each=ncol(factors[,1:4])),
  risk_premia = c(ifrp$risk_premia, aifrp$risk_premia, krs_frp$risk_premia),
  standard_errors = c(
    ifrp$standard_errors, aifrp$standard_errors, krs_frp$standard_errors
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
