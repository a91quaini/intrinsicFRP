# import package data on 15 risk factors and 42 test asset excess returns
factors = factors[,-1]
returns = returns[,-1]

# compute intrinsic factor risk premia and their standard errors
ifrp = IFRP(returns, factors, include_standard_errors = TRUE)

# set penalty parameters
penalty_parameters = seq(0., 1e-3, length.out = 1000)
# penalty_parameters = c(seq(0., 1e-4, length.out = 100), seq(1e-4, 1e-3, length.out = 100))

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
krs_frp = FRP(returns, factors)

PlotAdaptiveIFRPModelScore(aifrp, penalty_parameters)

