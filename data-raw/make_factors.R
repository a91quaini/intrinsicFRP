factors_ff5 <- read.csv("data-raw/F-F_Research_Data_5_Factors_2x3.csv")
factor_mom <- read.csv("data-raw/F-F_Momentum_Factor.csv")

start = factors_ff5[1,"Date"]
factor_mom = factor_mom[factor_mom[,"Date"] >= start,]

factors_ff5[,-1] = factors_ff5[,-1] / 100
factor_mom[,-1] = factor_mom[,-1] / 100

factors = cbind(
  factors_ff5[,-7],
  factor_mom[,-1]
)
factors = as.matrix(factors)
colnames(factors)[7] = colnames(factor_mom)[2]

usethis::use_data(factors, overwrite = TRUE)
