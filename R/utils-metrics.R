#' Mean squared error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar MSE.
#' @keywords internal
#' @noRd
.mse <- function(actual, predicted) {
  mean((actual - predicted)^2)
}

#' Root mean squared error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar RMSE.
#' @keywords internal
#' @noRd
.rmse <- function(actual, predicted) {
  sqrt(.mse(actual, predicted))
}

#' Mean absolute error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar MAE.
#' @keywords internal
#' @noRd
.mae <- function(actual, predicted) {
  mean(abs(actual - predicted))
}

#' Relative squared error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar RSE (SSE / SST).
#' @keywords internal
#' @noRd
.rse <- function(actual, predicted) {
  sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}
