# Global performance metrics: scores computed on the whole estimated
# composition vector of one bulk sample (one row per sample and algorithm),
# reconstitution scores when no ground truth is available, and the simplex
# KKT residual of the optimiser. Per-cell-type metrics live in
# 04_04_performance_metrics_celltype.R.

#' One row of composition-level scores for one bulk sample
#'
#' With a ground truth, compares \eqn{\hat{\boldsymbol{p}}} with
#' \eqn{\boldsymbol{p}^{\star}}; without it, compares the reconstituted bulk
#' \eqn{\boldsymbol{\mu}\hat{\boldsymbol{p}}} with \eqn{\boldsymbol{y}}.
#' Scores that do not apply are `NA`.
#'
#' @keywords internal
#' @noRd
.regression_global_row <- function(
  p_true,
  p_hat,
  y,
  mean_signature_matrix,
  sample_id,
  algorithm
) {
  if (is.null(p_true)) {
    y_hat <- as.vector(mean_signature_matrix %*% p_hat)
    return(tibble::tibble(
      sample_id = sample_id,
      algorithm = algorithm,
      tv = NA_real_,
      rmse = .rmse(y, y_hat),
      angular = NA_real_,
      sdid = NA_real_,
      maxae = NA_real_,
      reconstitution_mae = .mae(y, y_hat),
      reconstitution_cor = .pearson_safe(y, y_hat)
    ))
  }
  tibble::tibble(
    sample_id = sample_id,
    algorithm = algorithm,
    tv = .tv(p_true, p_hat),
    rmse = .rmse(p_true, p_hat),
    angular = .angular_distance(p_true, p_hat),
    sdid = .sdid(p_true, p_hat),
    maxae = .maxae(p_true, p_hat),
    reconstitution_mae = NA_real_,
    reconstitution_cor = NA_real_
  )
}

#' Root mean squared error
#'
#' @param actual Numeric vector of observed values.
#' @param predicted Numeric vector of predicted values.
#' @return Scalar RMSE.
#' @keywords internal
#' @noRd
.rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2))
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

#' Total variation / normalised MAE on the simplex
#'
#' @keywords internal
#' @noRd
.tv <- function(p, p_hat) {
  0.5 * sum(abs(as.numeric(p) - as.numeric(p_hat)))
}

#' \eqn{L_\infty} absolute error
#'
#' @keywords internal
#' @noRd
.maxae <- function(p, p_hat) {
  max(abs(as.numeric(p) - as.numeric(p_hat)))
}

#' Normalised angular distance on \eqn{[0, 1]}
#'
#' @keywords internal
#' @noRd
.angular_distance <- function(p, p_hat) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  denom <- sqrt(sum(p^2)) * sqrt(sum(p_hat^2))
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  cos_angle <- sum(p * p_hat) / denom
  cos_angle <- min(1, max(-1, cos_angle))
  (2 / pi) * acos(cos_angle)
}

#' HADACA-style SDID (\eqn{\sqrt{\sin\theta}})
#'
#' @keywords internal
#' @noRd
.sdid <- function(p, p_hat) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  denom <- sqrt(sum(p^2)) * sqrt(sum(p_hat^2))
  if (!is.finite(denom) || denom <= 0) {
    return(NA_real_)
  }
  cos_angle <- sum(p * p_hat) / denom
  sin_angle <- sqrt(1 - min(1, cos_angle^2))
  sqrt(sin_angle)
}

#' Safe Pearson correlation (NA if either vector is constant)
#'
#' @keywords internal
#' @noRd
.pearson_safe <- function(x, y) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 2L) {
    return(NA_real_)
  }
  sx <- stats::sd(x)
  sy <- stats::sd(y)
  if (!is.finite(sx) || !is.finite(sy) || sx <= 0 || sy <= 0) {
    return(NA_real_)
  }
  stats::cor(x, y, method = "pearson")
}

#' Simplex KKT / projected-score residual
#'
#' \eqn{\|\Pi_\Delta(p + \alpha g) - p\|_2 / \alpha} for maximisation,
#' with \eqn{\Pi_\Delta} the Euclidean projection onto the probability
#' simplex (sort-and-threshold algorithm).
#'
#' @keywords internal
#' @noRd
.kkt_residual <- function(p, grad, alpha = 1) {
  p <- as.numeric(p)
  grad <- as.numeric(grad)
  if (length(p) != length(grad) || length(p) == 0L) {
    return(NA_real_)
  }
  if (!is.finite(alpha) || alpha <= 0) {
    alpha <- 1
  }
  v <- p + alpha * grad
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u > (cssv - 1) / seq_along(u)))
  projected <- pmax(v - (cssv[[rho]] - 1) / rho, 0)
  sqrt(sum((projected - p)^2)) / alpha
}
