# Global performance metrics: scores computed on the whole estimated
# composition vector of one bulk sample (one row per sample and algorithm),
# reconstitution scores when no ground truth is available, and the simplex
# KKT residual of the optimiser. Per-cell-type metrics live in
# 04_04_performance_metrics_celltype.R.

#' One row of composition-level scores for one bulk sample
#'
#' With a ground truth, compares \eqn{\hat{\boldsymbol{p}}} with
#' \eqn{\boldsymbol{p}^{\star}}: total variation
#' \eqn{\tfrac12\lVert\hat p-p^{\star}\rVert_1}, RMSE, the largest absolute
#' error, the normalised angle \eqn{(2/\pi)\theta} and the HADACA score
#' SDID \eqn{\sqrt{\sin\theta}}, with \eqn{\theta} the angle between the
#' two compositions. Without a ground truth, compares the reconstituted
#' bulk \eqn{\boldsymbol{\mu}\hat{\boldsymbol{p}}} with \eqn{\boldsymbol{y}}
#' (RMSE, mean absolute error, Pearson correlation). Scores that do not
#' apply are `NA`.
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
    error <- y - y_hat
    return(tibble::tibble(
      sample_id = sample_id,
      algorithm = algorithm,
      tv = NA_real_,
      rmse = sqrt(mean(error^2)),
      angular = NA_real_,
      sdid = NA_real_,
      maxae = NA_real_,
      reconstitution_mae = mean(abs(error)),
      reconstitution_cor = .pearson_safe(y, y_hat)
    ))
  }
  p_true <- as.numeric(p_true)
  p_hat <- as.numeric(p_hat)
  error <- p_true - p_hat
  norms <- sqrt(sum(p_true^2) * sum(p_hat^2))
  cos_angle <- if (is.finite(norms) && norms > 0) {
    min(1, max(-1, sum(p_true * p_hat) / norms))
  } else {
    NA_real_
  }
  tibble::tibble(
    sample_id = sample_id,
    algorithm = algorithm,
    tv = 0.5 * sum(abs(error)),
    rmse = sqrt(mean(error^2)),
    angular = (2 / pi) * acos(cos_angle),
    sdid = sqrt(sqrt(1 - cos_angle^2)),
    maxae = max(abs(error)),
    reconstitution_mae = NA_real_,
    reconstitution_cor = NA_real_
  )
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
