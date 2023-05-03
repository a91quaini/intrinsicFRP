# Author: Alberto Quaini


#' Plot adaptive IFRP model score
#'
#' @name PlotAdaptiveIFRPModelScore
#' @description Plots the model score of the adaptive IFRP for each penalty
#' parameter value, highlighting the minimum attained model score and the
#' optimal one (which differs from the minimum in case the optimal parameter
#' was computed with the `one_stddev_rule` option TRUE).
#'
#' @param aifrp list containing the output of the function `OptimalAdaptiveIFRP`.
#' @param penalty_parameters n_parameters-dimensional vector of penalty
#' parameter values from smallest to largest.
#' @param legend_pos character vector indicating the legend position. Must be
#' one of "bottomright", "bottom", "bottomleft", "left", "topleft", "top",
#' "topright", "right" and "center" . Default is "bottomright".
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' penalty_parameters = seq(0., 1., length.out = 100)
#'
#' # compute optimal adaptive intrinsic factor risk premia and their standard errors
#' aifrp = OptimalAdaptiveIFRP(
#' returns,
#' factors,
#' penalty_parameters,
#' include_standard_errors = TRUE
#' )
#'
#' PlotAdaptiveIFRPModelScore(aifrp, penalty_parameters)
#'
#' @export
PlotAdaptiveIFRPModelScore = function(
  aifrp,
  penalty_parameters,
  legend_pos = "bottomright"
) {

  optimal_par = aifrp$penalty_parameter
  optimal_score = aifrp$model_score[penalty_parameters == aifrp$penalty_parameter]
  plot(
    penalty_parameters,
    aifrp$model_score,
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
    x = aifrp$penalty_parameter,
    y = optimal_score,
    pch = 15,
    col = "red"
  )
  graphics::abline(
    v = aifrp$penalty_parameter,
    col = "red",
    lty = 2
  )
  graphics::legend(
    legend_pos,
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

  stopifnot("`returns` does not contain numeric values" = is.numeric(returns))
  stopifnot("`factors` does not contain numeric values" = is.numeric(returns))
  stopifnot("`returns` and `factors` must have the same number of rows" = nrow(returns) == nrow(factors))
  stopifnot("`returns` contains more assets (columns) than observations (rows)" = nrow(returns) > ncol(returns))
  stopifnot("`factors` contains more variables (columns) than observations (rows)" = nrow(factors) > ncol(factors))
  stopifnot("the number of `returns` must be grater than the number of `factors`" = ncol(returns) > ncol(factors))
  stopifnot("`returns` must not contain missing values (NA/NaN)" = !anyNA(returns))
  stopifnot("`factors` must not contain missing values (NA/NaN)" = !anyNA(factors))

}
