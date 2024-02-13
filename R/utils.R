# Author: Alberto Quaini


# Plot Oracle TFRP model score
#
# @name PlotOracleTFRPModelScore
# @description Plots the model score of the Oracle tradable factor risk premia
# for each penalty
# parameter value, highlighting the minimum attained model score and the
# optimal one (which differs from the minimum in case the optimal parameter
# was computed with the `one_stddev_rule` option TRUE).
#
# @param oracle_tfrp list containing the output of the function `OracleTFRP`.
# @param penalty_parameters n_parameters-dimensional vector of penalty
# parameter values from smallest to largest.
# @param legend_pos character vector indicating the legend position. Must be
# one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top",
# "topright", "right" and "center" . Default is "bottomright".
#
PlotOracleTFRPModelScore = function(
  oracle_tfrp,
  penalty_parameters
) {

  optimal_score = oracle_tfrp$model_score[penalty_parameters == oracle_tfrp$penalty_parameter]

  plot(
    penalty_parameters,
    oracle_tfrp$model_score,
    frame = FALSE,
    type = "b",
    pch = 20,
    lty = 2,
    #col = "black",
    xlab = "Penalty parameter",
    ylab = "Model score",
    cex = 0.6,
    cex.lab = 1.3
  )

  graphics::points(
    x = oracle_tfrp$penalty_parameter,
    y = optimal_score,
    pch = 15,
    col = "red"
  )

  graphics::abline(
    v = oracle_tfrp$penalty_parameter,
    col = "red",
    lty = 2
  )

  graphics::legend(
    "topright",
    legend = c("Optimal"),
    col = "red",
    pch = 15,
    cex = 1.1
  )

}

## Check returns and factors
# Checks that the main data arguments, namely returns and factors, are
# conforming to the packgae implementation.
CheckData = function(returns, factors) {

  stopifnot("`returns` must contain numeric values" = is.numeric(returns))
  stopifnot("`factors` must contain numeric values" = is.numeric(factors))
  stopifnot("`returns` and `factors` must have the same number of rows" = nrow(returns) == nrow(factors))
  stopifnot("`returns` contains more assets (columns) than observations (rows)" = nrow(returns) > ncol(returns))
  stopifnot("`factors` contains more variables (columns) than observations (rows)" = nrow(factors) > ncol(factors))
  stopifnot("the number of `returns` must be grater or equal to the number of `factors`" = ncol(returns) >= ncol(factors))
  stopifnot("`returns` must not contain missing values (NA/NaN)" = !anyNA(returns))
  stopifnot("`factors` must not contain missing values (NA/NaN)" = !anyNA(factors))

}
