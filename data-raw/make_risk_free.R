factors_ff5 <- read.csv("data-raw/F-F_Research_Data_5_Factors_2x3.csv")

risk_free = as.matrix(factors_ff5[,c(1, 7)])
risk_free[,2] = risk_free[,2] / 100

usethis::use_data(risk_free, overwrite = TRUE)
