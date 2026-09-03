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
  sst <- sum((actual - mean(actual))^2)
  if (!is.finite(sst) || sst <= 0) {
    return(NA_real_)
  }
  sum((actual - predicted)^2) / sst
}

#' Euclidean projection onto the probability simplex
#'
#' @param v Numeric vector.
#' @return Non-negative vector summing to 1.
#' @keywords internal
#' @noRd
.project_simplex <- function(v) {
  v <- as.numeric(v)
  n <- length(v)
  if (n == 0L) {
    stop("`v` must be non-empty.", call. = FALSE)
  }
  u <- sort(v, decreasing = TRUE)
  cssv <- cumsum(u)
  rho <- max(which(u > (cssv - 1) / seq_len(n)))
  theta <- (cssv[[rho]] - 1) / rho
  pmax(v - theta, 0)
}

#' Simplex KKT / projected-score residual
#'
#' \eqn{\|\Pi_\Delta(p + \alpha g) - p\|_2 / \alpha} for maximisation.
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
  step <- .project_simplex(p + alpha * grad) - p
  sqrt(sum(step^2)) / alpha
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
  if (length(x) < 2L || length(y) < 2L) {
    return(NA_real_)
  }
  if (stats::sd(x) <= 0 || stats::sd(y) <= 0) {
    return(NA_real_)
  }
  stats::cor(x, y, method = "pearson")
}

#' Presence / spillover diagnostics for one composition pair
#'
#' @keywords internal
#' @noRd
.presence_counts <- function(p, p_hat, threshold = 1e-4) {
  p <- as.numeric(p)
  p_hat <- as.numeric(p_hat)
  present_true <- p > threshold
  present_hat <- p_hat > threshold
  list(
    tp = sum(present_true & present_hat),
    fp = sum(!present_true & present_hat),
    fn = sum(present_true & !present_hat),
    tn = sum(!present_true & !present_hat),
    false_positive_mass = sum(p_hat[!present_true])
  )
}

#' F1 from presence counts
#'
#' @keywords internal
#' @noRd
.f1_from_counts <- function(tp, fp, fn) {
  if ((tp + fp + fn) == 0L) {
    return(NA_real_)
  }
  precision <- if ((tp + fp) == 0L) NA_real_ else tp / (tp + fp)
  recall <- if ((tp + fn) == 0L) NA_real_ else tp / (tp + fn)
  if (
    !is.finite(precision) || !is.finite(recall) || (precision + recall) == 0
  ) {
    return(0)
  }
  2 * precision * recall / (precision + recall)
}

#' Process memory in bytes (PSS on Linux when `ps` is available)
#'
#' @keywords internal
#' @noRd
.process_memory_bytes <- function() {
  if (!requireNamespace("ps", quietly = TRUE)) {
    return(NA_real_)
  }
  info <- tryCatch(
    ps::ps_memory_full_info(),
    error = function(e) {
      tryCatch(ps::ps_memory_info(), error = function(e2) NULL)
    }
  )
  if (is.null(info)) {
    return(NA_real_)
  }
  pss <- info[["pss"]]
  if (!is.null(pss) && is.finite(pss)) {
    return(as.numeric(pss))
  }
  rss <- info[["rss"]]
  if (!is.null(rss) && is.finite(rss)) {
    return(as.numeric(rss))
  }
  NA_real_
}

#' Monte Carlo SE of a mean
#'
#' @keywords internal
#' @noRd
.mcse_mean <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 2L) {
    return(NA_real_)
  }
  stats::sd(x) / sqrt(n)
}

#' Monte Carlo SE of an empirical coverage probability
#'
#' @keywords internal
#' @noRd
.mcse_coverage <- function(covered) {
  coverage_mc_interval(covered, method = "wald")$mcse
}

#' Monte Carlo interval for an empirical coverage probability
#'
#' Interval for a binomial coverage rate
#' \eqn{\hat\pi=X/N} from independent coverage indicators.
#' Wilson is the default because the parameter is bounded in
#' \eqn{[0,1]} and Wald intervals behave poorly near 0 or 1.
#'
#' @param covered Logical or 0/1 coverage indicators (one per replicate).
#' @param conf_level Confidence level for the interval around
#'   \eqn{\hat\pi} (default `0.95`).
#' @param method `"wilson"` (default), `"wald"`, or `"agresti_coull"`.
#'
#' @return A list with `n`, `successes`, `coverage`, `mcse`, `lower`,
#'   `upper`, and `method`.
#'
#' @examples
#' coverage_mc_interval(c(TRUE, TRUE, TRUE, FALSE))
#' @export
coverage_mc_interval <- function(
  covered,
  conf_level = 0.95,
  method = c("wilson", "wald", "agresti_coull")
) {
  method <- .match_arg_case_insensitive(
    method,
    c(
      "wilson",
      "wald",
      "agresti_coull"
    )
  )
  covered <- as.numeric(covered)
  covered <- covered[is.finite(covered)]
  n <- length(covered)
  if (n == 0L) {
    return(list(
      n = 0L,
      successes = 0L,
      coverage = NA_real_,
      mcse = NA_real_,
      lower = NA_real_,
      upper = NA_real_,
      method = method
    ))
  }
  x <- sum(covered)
  p_hat <- x / n
  mcse <- sqrt(p_hat * (1 - p_hat) / n)
  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  bounds <- switch(
    method,
    wald = {
      c(
        max(0, p_hat - z * mcse),
        min(1, p_hat + z * mcse)
      )
    },
    wilson = {
      d <- 1 + z^2 / n
      centre <- (p_hat + z^2 / (2 * n)) / d
      half <- z /
        d *
        sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
      c(max(0, centre - half), min(1, centre + half))
    },
    agresti_coull = {
      n_tilde <- n + z^2
      p_tilde <- (x + z^2 / 2) / n_tilde
      half <- z * sqrt(p_tilde * (1 - p_tilde) / n_tilde)
      c(max(0, p_tilde - half), min(1, p_tilde + half))
    }
  )
  list(
    n = n,
    successes = as.integer(x),
    coverage = p_hat,
    mcse = mcse,
    lower = bounds[[1L]],
    upper = bounds[[2L]],
    method = method
  )
}

#' Empty three-block benchmark list
#'
#' @keywords internal
#' @noRd
.empty_benchmark_metrics <- function() {
  list(
    regression = list(
      global = tibble::tibble(),
      cell_type = tibble::tibble()
    ),
    monte_carlo = tibble::tibble(),
    optimisation = tibble::tibble()
  )
}

#' Align a proportion vector with signature colnames
#'
#' @keywords internal
#' @noRd
.as_named_p <- function(p, nms) {
  p <- as.numeric(p)
  if (!is.null(nms) && length(nms) == length(p)) {
    names(p) <- nms
  }
  p
}

#' Coerce a solver return value to a named proportion vector
#'
#' @keywords internal
#' @noRd
.coerce_estimated_p <- function(out, nms) {
  if (is.list(out) && !is.null(out$coefficients)) {
    p <- out$coefficients
  } else {
    p <- out
  }
  if (is.matrix(p) || is.data.frame(p)) {
    p <- as.numeric(p[, 1L, drop = TRUE])
  }
  .as_named_p(p, nms)
}

#' Optional standard errors from a solver return value
#'
#' @keywords internal
#' @noRd
.extract_estimated_se <- function(out, n_celltypes) {
  if (!is.list(out)) {
    return(rep(NA_real_, n_celltypes))
  }
  if (!is.null(out$se)) {
    return(as.numeric(out$se)[seq_len(n_celltypes)])
  }
  if (!is.null(out$vcov)) {
    return(sqrt(pmax(diag(as.matrix(out$vcov)), 0)))
  }
  rep(NA_real_, n_celltypes)
}

#' Map over bulk samples, optionally with `furrr`
#'
#' @keywords internal
#' @noRd
.map_samples <- function(n, fun, cores) {
  idx <- seq_len(n)
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    stop("`cores` must be a positive integer.", call. = FALSE)
  }
  if (cores == 1L) {
    return(lapply(idx, fun))
  }
  .check_suggested_package("furrr", "deconvolute_ratios")
  .check_suggested_package("future", "deconvolute_ratios")
  old_plan <- future::plan(future::multisession, workers = cores)
  on.exit(future::plan(old_plan), add = TRUE)
  furrr::future_map(
    idx,
    fun,
    .options = furrr::furrr_options(
      # L'Ecuyer-CMRG streams per worker (R >= 4.1 / future.seed).
      seed = TRUE,
      packages = "DeCovarT"
    )
  )
}
