# Benchmark support for the two likelihoods and their confidence intervals:
# the per-sample diagnostics of 04_01_run_benchmark.R need the log-likelihood,
# score and curvature of the model that generated the bulk, and the interval
# table compares the interval menu of 03_01 / 03_04 on the same fits.

#' Likelihood functions of the benchmarked model
#'
#' `population` is the convolution \eqn{\mathcal{N}(\boldsymbol{\mu}
#' \boldsymbol{p},\sum_j p_j^2\boldsymbol{\Sigma}_j)}; `single_cell` is the
#' mixture of cells with its two scales profiled out, so every function is
#' evaluated at the closed-form scales of the composition.
#'
#' These are the functions [run_simulation_benchmark()] uses to score an
#' estimate (log-likelihood at the estimate and at the truth, KKT residual
#' through [kkt_residual()]); they let estimates computed outside the
#' package be scored on the same terms.
#'
#' @param model `"population"` or `"single_cell"`.
#' @param mean_signature_matrix \eqn{G\times J} mean matrix.
#' @param Sigma \eqn{G\times G\times J} covariance array.
#'
#' @return A list of functions of `(p, y)`: `loglik`, `score`, `se` (Wald
#'   standard errors of the composition) and `curvature` (local-maximum
#'   flag and eigenvalue range of the Hessian in ILR coordinates, `NA` on
#'   a face).
#'
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' model <- benchmark_likelihood("population", mu, Sigma)
#' y <- drop(mu %*% c(0.4, 0.6))
#' model$loglik(c(0.4, 0.6), y)
#' kkt_residual(c(0.4, 0.6), model$score(c(0.4, 0.6), y))
#' @export
benchmark_likelihood <- function(
  model = c("population", "single_cell"),
  mean_signature_matrix,
  Sigma
) {
  model <- match.arg(model)
  no_curvature <- list(
    local_maximum = NA,
    min_eigenvalue = NA_real_,
    max_eigenvalue = NA_real_
  )
  if (identical(model, "population")) {
    return(list(
      loglik = function(p, y) {
        loglik_multivariate(p, y, mean_signature_matrix, Sigma)
      },
      score = function(p, y) {
        gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma)
      },
      se = function(p, y) {
        .ilr_wald_se(p, mean_signature_matrix, Sigma, warn = FALSE)
      },
      curvature = function(p, y) {
        if (min(p) <= 0) {
          return(no_curvature)
        }
        diagnostics <- boundary_diagnostics(p, y, mean_signature_matrix, Sigma)
        as.list(diagnostics[c(
          "local_maximum",
          "min_eigenvalue",
          "max_eigenvalue"
        )])
      }
    ))
  }
  scales_at <- function(p, y) {
    .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p, strict = FALSE)
  }
  list(
    loglik = function(p, y) {
      scales <- scales_at(p, y)
      if (is.null(scales)) {
        return(NA_real_)
      }
      loglik_single_cell(
        p,
        y,
        mean_signature_matrix,
        Sigma,
        scales$mean_scale,
        scales$variance_scale
      )
    },
    score = function(p, y) {
      scales <- scales_at(p, y)
      gradient_loglik_single_cell(
        p,
        y,
        mean_signature_matrix,
        Sigma,
        scales$mean_scale,
        scales$variance_scale
      )
    },
    se = function(p, y) {
      vcov <- vcov_single_cell(
        p,
        y,
        mean_signature_matrix,
        Sigma,
        warn = FALSE
      )$vcov
      stats::setNames(sqrt(pmax(diag(as.matrix(vcov)), 0)), names(p))
    },
    curvature = function(p, y) {
      if (min(p) <= 0) {
        return(no_curvature)
      }
      # Same criterion as boundary_diagnostics(): small ILR score and a
      # negative-definite Hessian of the profile log-likelihood.
      z <- isometric_log_ratio(p)
      score <- gradient_profile_loglik_single_cell_ilr(
        z,
        y,
        mean_signature_matrix,
        Sigma
      )
      ev <- eigen(
        hessian_profile_loglik_single_cell_ilr(
          z,
          y,
          mean_signature_matrix,
          Sigma
        ),
        symmetric = TRUE,
        only.values = TRUE
      )$values
      list(
        local_maximum = sqrt(sum(score^2)) < 1e-4 && max(ev) < 0,
        min_eigenvalue = min(ev),
        max_eigenvalue = max(ev)
      )
    }
  )
}

#' Confidence intervals of one fitted bulk sample
#'
#' Evaluates the interval menu at the estimate `p_hat` of the benchmarked
#' algorithm:
#' * `wald`: delta method on \eqn{\boldsymbol{p}}, symmetric;
#' * `ilr`: Wald interval in ILR coordinates mapped back to the simplex;
#' * `profile`: likelihood-ratio inversion, chart-invariant (it refits);
#' * `bootstrap`: percentile interval of a parametric bootstrap refitted by
#'   `ci$optimiser`;
#' * `oracle`: `p_hat` plus or minus the Wald half-width evaluated at the
#'   true composition, the width an analyst would get with the true
#'   information.
#'
#' For the single-cell model, `_reml` appended to `wald`, `ilr` or
#' `profile` applies the REML calibration (\eqn{\nu=G-J}). Single-cell
#' Wald, ILR and bootstrap intervals also cover the scales \eqn{c} and
#' \eqn{\phi} (rows `mean_scale` and `variance_scale`) and the derived
#' \eqn{a=\phi/c} and \eqn{n=c^{2}/\phi} (rows `platform_scale_a` and
#' `cell_count_n`). A method that
#' fails, or is undefined because `p_hat` lies on a face, gives `NA`
#' bounds.
#'
#' @param ci List with `methods`, `level`, `n_boot`, `n_draws` and
#'   `optimiser`.
#' @return A named list of two-column matrices (lower, upper), one per
#'   method, with one row per parameter.
#' @noRd
.benchmark_intervals <- function(
  y,
  p_hat,
  p_true,
  mean_signature_matrix,
  Sigma,
  ci,
  model
) {
  level <- if (is.null(ci$level)) 0.95 else ci$level
  n_boot <- if (is.null(ci$n_boot)) 200L else ci$n_boot
  n_draws <- if (is.null(ci$n_draws)) 2000L else ci$n_draws
  optimiser <- if (is.null(ci$optimiser)) "newton" else ci$optimiser
  nms <- names(p_hat)
  bounds <- function(x) {
    x <- as.matrix(x)[, seq_len(2L), drop = FALSE]
    dimnames(x) <- list(rownames(x), c("lower", "upper"))
    x
  }
  missing_bounds <- bounds(matrix(
    NA_real_,
    length(nms),
    2L,
    dimnames = list(nms, NULL)
  ))
  population <- identical(model, "population")
  one_method <- function(method) {
    reml <- grepl("_reml$", method)
    base <- sub("_reml$", "", method)
    calibration <- .single_cell_calibration(
      if (reml) "reml" else "ml",
      n_genes = nrow(mean_signature_matrix),
      n_types = length(nms)
    )
    if (base %in% c("wald", "ilr", "oracle") && min(p_hat) <= 0) {
      return(missing_bounds)
    }
    if (population) {
      switch(
        base,
        wald = .wald_bounds(
          p_hat,
          .ilr_wald_se(p_hat, mean_signature_matrix, Sigma),
          level
        ),
        oracle = .wald_bounds(
          p_hat,
          .ilr_wald_se(p_true, mean_signature_matrix, Sigma),
          level
        ),
        ilr = suppressWarnings(confint_ilr_backtransform(
          y,
          mean_signature_matrix,
          Sigma,
          p = p_hat,
          level = level,
          n_draws = n_draws
        ))$interval,
        profile = confint_profile_decovart(
          y,
          mean_signature_matrix,
          Sigma,
          level = level
        )[, c("lower", "upper"), drop = FALSE],
        bootstrap = bootstrap_decovart(
          y,
          mean_signature_matrix,
          Sigma,
          p = p_hat,
          n_boot = n_boot,
          level = level,
          optimiser = optimiser
        )$interval,
        stop("Unknown interval method `", method, "`.", call. = FALSE)
      )
    } else {
      switch(
        base,
        wald = ,
        ilr = {
          pieces <- .composition_intervals(
            p_hat,
            y,
            mean_signature_matrix,
            Sigma,
            level,
            how = base,
            n_draws = n_draws,
            seed = NULL,
            calibration = calibration
          )
          rbind(bounds(pieces$interval), bounds(pieces$nuisance$interval))
        },
        oracle = {
          truth <- vcov_single_cell(
            p_true,
            y,
            mean_signature_matrix,
            Sigma,
            warn = FALSE
          )
          .wald_bounds(p_hat, sqrt(pmax(diag(truth$vcov), 0)), level)
        },
        profile = confint_profile_single_cell(
          y,
          mean_signature_matrix,
          Sigma,
          level = level,
          variance = calibration$variance
        )[, c("lower", "upper"), drop = FALSE],
        bootstrap = {
          boot <- bootstrap_single_cell(
            y,
            mean_signature_matrix,
            Sigma,
            n_boot = n_boot,
            level = level,
            optimiser = optimiser
          )
          rbind(bounds(boot$interval), bounds(boot$scale_interval))
        },
        stop("Unknown interval method `", method, "`.", call. = FALSE)
      )
    }
  }
  stats::setNames(
    lapply(ci$methods, function(method) {
      out <- tryCatch(
        suppressWarnings(one_method(method)),
        error = function(e) missing_bounds
      )
      bounds(out)
    }),
    ci$methods
  )
}

#' Coverage and width of every interval method over the Monte Carlo samples
#'
#' @param intervals List over samples of the output of
#'   [.benchmark_intervals()].
#' @param p_hat,p_true \eqn{J\times N} estimates and targets.
#' @param nuisance_hat,nuisance_true Optional \eqn{2\times N} matrices of
#'   the estimated and true scales (`mean_scale`, `variance_scale`); the
#'   platform scale \eqn{a=\phi/c} and the cell count \eqn{n=c^{2}/\phi}
#'   are derived from them and scored as well.
#' @return A tibble with one row per algorithm, interval method and
#'   parameter.
#' @noRd
.interval_table <- function(
  intervals,
  p_hat,
  p_true,
  nuisance_hat = NULL,
  nuisance_true = NULL,
  coverage_interval = "wilson",
  algorithm = NA_character_
) {
  if (length(intervals) == 0L || all(vapply(intervals, is.null, TRUE))) {
    return(tibble::tibble())
  }
  derive <- function(scales) {
    if (is.null(scales)) NULL else .derived_scales(scales)
  }
  estimate <- rbind(p_hat, derive(nuisance_hat))
  truth <- rbind(p_true, derive(nuisance_true))
  methods <- unique(unlist(lapply(intervals, names)))
  composition <- rownames(p_hat)
  purrr::map_dfr(methods, function(method) {
    parameters <- intersect(
      rownames(truth),
      unique(unlist(lapply(intervals, function(x) rownames(x[[method]]))))
    )
    purrr::map_dfr(parameters, function(parameter) {
      bound <- function(side) {
        vapply(
          intervals,
          function(x) {
            # A bulk whose interval failed may lack this parameter's row.
            m <- x[[method]]
            if (
              is.null(m) ||
                !parameter %in% rownames(m) ||
                !side %in% colnames(m)
            ) {
              return(NA_real_)
            }
            value <- m[parameter, side]
            if (length(value) != 1L) NA_real_ else value
          },
          numeric(1)
        )
      }
      lower <- bound("lower")
      upper <- bound("upper")
      target <- truth[parameter, ]
      finite <- is.finite(lower) & is.finite(upper)
      width <- upper - lower
      rate <- coverage_mc_interval(
        (lower <= target & target <= upper)[finite],
        method = coverage_interval
      )
      tibble::tibble(
        algorithm = algorithm,
        ci_method = method,
        parameter = parameter,
        is_composition = parameter %in% composition,
        n_intervals = sum(finite),
        coverage = rate$coverage,
        coverage_lower = rate$lower,
        coverage_upper = rate$upper,
        mean_width = if (any(finite)) mean(width[finite]) else NA_real_,
        median_width = if (any(finite)) {
          stats::median(width[finite])
        } else {
          NA_real_
        },
        mean_lower = if (any(finite)) mean(lower[finite]) else NA_real_,
        mean_upper = if (any(finite)) mean(upper[finite]) else NA_real_,
        mean_estimate = mean(estimate[parameter, ], na.rm = TRUE),
        mean_truth = mean(target, na.rm = TRUE),
        outside_unit = if (parameter %in% composition && any(finite)) {
          mean((lower < 0 | upper > 1)[finite])
        } else {
          NA_real_
        }
      )
    })
  })
}
