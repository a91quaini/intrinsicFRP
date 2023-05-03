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
#'
#' @examples
#' # import package data on 15 risk factors and 42 test asset excess returns
#' factors = intrinsicFRP::factors[,-1]
#' returns = intrinsicFRP::returns[,-1]
#'
#' penalty_parameters = seq(0., 1., length.out = 100)
#'
#' # compute optimal adaptive intrinsic factor risk premia and their standard errors
# #' aifrp = OptimalAdaptiveIFRP(
# #' returns,
# #' factors,
# #' penalty_parameters,
# #' include_standard_errors = TRUE
# #' )
# #'
# #' PlotAdaptiveIFRPModelScore(aifrp, penalty_parameters)
#'
#' @export
PlotAdaptiveIFRPModelScore = function(
  aifrp,
  penalty_parameters
) {

  idx_min_score = which.min(aifrp$model_score)
  sd_score = stats::sd(aifrp$model_score)
  idx_optimal = which(penalty_parameters == aifrp$penalty_parameter)

  df = data.frame(
    model_score = aifrp$model_score,
    penalty_parameters = penalty_parameters
  )

  ggplot2::ggplot(df, ggplot2::aes(x = .data$penalty_parameters, y = .data$model_score)) +
    # ggplot2::theme_bw() +
    ggplot2::theme(
      text=ggplot2::element_text(size=22),
      panel.background = ggplot2::element_rect(fill = "white"),
      panel.grid.major.y = ggplot2::element_line(
        color="gray", linewidth =0.25
      )
    ) +
    ggplot2::xlab("Penalty parameter") +
    ggplot2::ylab("Model score") +
    ggplot2::geom_line() +
    ggplot2::geom_point(
      ggplot2::aes(
        x = .data$penalty_parameters[idx_min_score],
        y = .data$model_score[idx_min_score]
      ),
      color="blue",
      size = 2.5
    ) +
    ggplot2::geom_label(
      label = "Minimum",
      x = penalty_parameters[idx_min_score],
      y = aifrp$model_score[idx_min_score],
      hjust = 0, vjust = -1,
      color = "blue"
    ) +
    ggplot2::geom_hline(
      yintercept = aifrp$model_score[idx_min_score],
      linetype = "dashed",
      color = "blue"
    ) +
    ggplot2::geom_point(ggplot2::aes(
      x = aifrp$penalty_parameter,
      y = .data$model_score[idx_optimal]),
      color="purple",
      size = 2.5
    ) +
    ggplot2::geom_label(
      label = "Optimal",
      x = aifrp$penalty_parameter,
      y = aifrp$model_score[idx_optimal],
      hjust = 0, vjust = 2,
      color = "purple"
    ) +
    ggplot2::geom_hline(
      yintercept = aifrp$model_score[idx_optimal],
      linetype = "dashed",
      color = "purple"
    )


}

# PlotAdaptiveIFRPModelScore1 = function(aifrp, penalty_parameters) {
#
#   plot(penalty_parameters, aifrp$model_score)
#
# }

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
