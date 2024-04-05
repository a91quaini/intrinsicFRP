## obtain data on portfolio excess returns from the Kenneth French data library
# from csv file to data.frame
# assumes that superfluous text and equally weighted returns were removed from
# the csv files, and the first column was named "Date"
# also, requires data on FF3 factors (for the risk-free return)

path = "data-raw/"
ff5  = "F-F_Research_Data_5_Factors_2x3.csv"
mebeme = "25_Portfolios_5x5.CSV"
ind = "17_Industry_Portfolios.CSV"

# source data on the risk free and portfolio returns
factors_ff5 = read.csv(file=paste(path, ff5, sep = ""))
returns_mebeme = read.csv(paste(path, mebeme, sep = ""))
returns_ind = read.csv(paste(path, ind, sep = ""))

start = factors_ff5[1,"Date"]
returns_mebeme = returns_mebeme[returns_mebeme[,"Date"] >= start,]
returns_ind = returns_ind[returns_ind[,"Date"] >= start,]
returns = cbind(
  returns_mebeme,
  returns_ind[,-1]
)
RF = factors_ff5[,7]

# subtract the risk free and divide by 100 (to obtain excess returns not in percentage)
returns[, -1] = (returns[, -1] - RF) / 100
returns = as.matrix(returns)

usethis::use_data(returns, overwrite = TRUE)
