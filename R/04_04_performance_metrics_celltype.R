# Cell-type performance metrics: scores computed for each cellular ratio
# p_j across the bulk samples (Monte Carlo replicates) of one scenario and
# algorithm: Pearson correlation, presence F1 and false-positive mass, bias,
# standard deviations, RMSE, and the coverage rate with its Monte Carlo
# interval. Composition-level metrics live in
# 04_03_performance_metrics_global.R.

#' Across-sample Pearson correlation and presence scores per cell type
#'
#' A cell type is present when its proportion exceeds
#' `presence_threshold`; failed fits (`NA`) count as absent. F1 is
#' \eqn{2\mathrm{TP}/(2\mathrm{TP}+\mathrm{FP}+\mathrm{FN})}, `NA` when the
#' type is never present in truth or estimate.
#'
#' @keywords internal
#' @noRd
.regression_cell_type_table <- function(
  p_true,
  p_hat,
  presence_threshold,
  algorithm
) {
  if (is.null(p_true)) {
    return(tibble::tibble(
      algorithm = character(),
      cell_type = character(),
      pearson = numeric(),
      presence_f1 = numeric(),
      false_positive_mass = numeric()
    ))
  }
  purrr::map_dfr(seq_len(nrow(p_hat)), function(j) {
    truth <- p_true[j, ]
    estimate <- p_hat[j, ]
    truth[!is.finite(truth)] <- 0
    estimate[!is.finite(estimate)] <- 0
    present_true <- truth > presence_threshold
    present_hat <- estimate > presence_threshold
    tp <- sum(present_true & present_hat)
    errors <- sum(!present_true & present_hat) +
      sum(present_true & !present_hat)
    tibble::tibble(
      algorithm = algorithm,
      cell_type = rownames(p_hat)[[j]],
      pearson = .pearson_safe(p_true[j, ], p_hat[j, ]),
      presence_f1 = if (tp + errors == 0) {
        NA_real_
      } else {
        2 * tp / (2 * tp + errors)
      },
      false_positive_mass = sum(estimate[!present_true]) / ncol(p_hat)
    )
  })
}

#' @keywords internal
#' @noRd
.monte_carlo_table <- function(
  p_true,
  p_hat,
  se,
  lower,
  upper,
  coverage_interval,
  algorithm
) {
  cell_names <- rownames(p_hat)
  if (is.null(p_true)) {
    return(tibble::tibble(
      algorithm = character(),
      cell_type = character(),
      bias = numeric(),
      empirical_sd = numeric(),
      mean_model_sd = numeric(),
      mean_model_se = numeric(),
      se_sd_ratio = numeric(),
      theoretical_se = numeric(),
      rmse = numeric(),
      coverage = numeric(),
      coverage_lower = numeric(),
      coverage_upper = numeric(),
      coverage_interval = character(),
      mean_interval_width = numeric(),
      mcse_coverage = numeric()
    ))
  }
  n_celltypes <- nrow(p_hat)
  n_samples <- ncol(p_hat)
  purrr::map_dfr(seq_len(n_celltypes), function(j) {
    err <- p_hat[j, ] - p_true[j, ]
    se_j <- se[j, ]
    empirical_sd <- if (n_samples < 2L) {
      NA_real_
    } else {
      stats::sd(p_hat[j, ], na.rm = TRUE)
    }
    covered <- lower[j, ] <= p_true[j, ] & p_true[j, ] <= upper[j, ]
    interval <- coverage_mc_interval(
      covered,
      method = coverage_interval
    )
    se_finite <- se_j[is.finite(se_j)]
    mean_model_se <- if (length(se_finite) == 0L) {
      NA_real_
    } else {
      mean(se_finite)
    }
    mean_model_sd <- if (length(se_finite) == 0L) {
      NA_real_
    } else {
      sqrt(mean(se_finite^2))
    }
    width <- mean(upper[j, ] - lower[j, ], na.rm = TRUE)
    if (!is.finite(width)) {
      width <- NA_real_
    }
    tibble::tibble(
      algorithm = algorithm,
      cell_type = cell_names[[j]],
      bias = mean(err, na.rm = TRUE),
      empirical_sd = empirical_sd,
      mean_model_sd = mean_model_sd,
      mean_model_se = mean_model_se,
      se_sd_ratio = mean_model_se / empirical_sd,
      theoretical_se = NA_real_,
      rmse = sqrt(mean(err^2, na.rm = TRUE)),
      coverage = interval$coverage,
      coverage_lower = interval$lower,
      coverage_upper = interval$upper,
      coverage_interval = interval$method,
      mean_interval_width = width,
      mcse_coverage = interval$mcse
    )
  })
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
