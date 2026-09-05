#' Validate generative-model parameters \eqn{\theta}
#'
#' Checks that `true_theta` stores the DeCovarT mixture parameters:
#' * `p`: mixing proportions of length \eqn{J}, or a matrix
#'   \eqn{J\times N} of sample-wise proportions (each column on the simplex);
#' * `mu`: mean signature \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}
#'   (also accepts \eqn{J\times G});
#' * `sigma` and/or `Theta`: array
#'   \eqn{G\times G\times J} of covariances and/or precision matrices.
#'
#' @param true_theta Named list with the fields above.
#' @param require_p Logical; if `TRUE` (default), `p` must be present and
#'   valid. If `FALSE`, a missing `p` is allowed (callers may default it).
#' @param J Optional expected number of cell types. When `NULL`, inferred
#'   from `mu` / the third dimension of `sigma` or `Theta`.
#' @param second_moment Which second-moment field to require:
#'   `"sigma"`, `"Theta"`, or `"either"` (default: at least one of
#'   `sigma` / `Theta`). Matching is case-insensitive.
#'
#' @srrstats {G2.3} Restricted character input (`second_moment`).
#' @srrstats {G2.3a} Validated via `.match_arg_case_insensitive()` (a
#'   `match.arg()` equivalent).
#' @srrstats {G2.3b} Matching is case-insensitive (`tolower()`).
#'
#' @return `TRUE` invisibly if the structure is valid; otherwise stops
#'   with an informative error.
#'
#' @seealso [compute_average_overlap()], [compute_average_jeffreys()]
#' @export
#' @examples
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' check_true_theta(theta)
check_true_theta <- function(
  true_theta,
  require_p = TRUE,
  J = NULL,
  second_moment = c("either", "sigma", "Theta")
) {
  .parse_true_theta(
    true_theta = true_theta,
    require_p = require_p,
    J = J,
    second_moment = second_moment
  )
  invisible(TRUE)
}

#' Parse and standardise \eqn{\theta} (internal)
#'
#' @return Named list with `p` (length-\eqn{J} vector or `NULL`),
#'   `mu` (\eqn{G\times J}), `sigma` and/or `Theta` (\eqn{G\times G\times J}),
#'   `G`, and `J`.
#' @keywords internal
#' @noRd
.parse_true_theta <- function(
  true_theta,
  require_p = TRUE,
  J = NULL,
  second_moment = c("either", "sigma", "Theta")
) {
  second_moment <- .match_arg_case_insensitive(
    second_moment,
    c("either", "sigma", "Theta")
  )
  if (!is.list(true_theta)) {
    stop("`true_theta` must be a list.", call. = FALSE)
  }
  if (is.null(true_theta$mu)) {
    stop("`true_theta` must contain `mu`.", call. = FALSE)
  }

  # --- Second moments: require Sigma and/or precision Theta -------------------
  # Each must be a G x G x J array (one PD matrix per cell type). When both are
  # supplied they must share the same (G, J).
  has_sigma <- !is.null(true_theta$sigma)
  has_theta <- !is.null(true_theta$Theta)
  if (identical(second_moment, "sigma") && !has_sigma) {
    stop("`true_theta` must contain `sigma`.", call. = FALSE)
  }
  if (identical(second_moment, "Theta") && !has_theta) {
    stop("`true_theta` must contain `Theta`.", call. = FALSE)
  }
  if (identical(second_moment, "either") && !has_sigma && !has_theta) {
    stop(
      "`true_theta` must contain `sigma` and/or `Theta`.",
      call. = FALSE
    )
  }

  # Helper: enforce cubic layout G x G x J and return (G, J).
  # Shared implementation lives in utils.R as `.assert_ggj_array()`.
  .check_ggj_array <- function(arr, nm) {
    .assert_ggj_array(arr, paste0("true_theta$", nm))
  }

  dims_list <- list()
  if (has_sigma) {
    dims_list$sigma <- .check_ggj_array(true_theta$sigma, "sigma")
  }
  if (has_theta) {
    dims_list$Theta <- .check_ggj_array(true_theta$Theta, "Theta")
  }
  if (length(dims_list) == 2L) {
    if (
      !identical(dims_list$sigma$G, dims_list$Theta$G) ||
        !identical(dims_list$sigma$J, dims_list$Theta$J)
    ) {
      stop("`sigma` and `Theta` must share the same G x G x J dims.")
    }
  }
  G <- dims_list[[1L]]$G
  n_celltypes <- dims_list[[1L]]$J

  # --- Mean signature mu: genes x cell types ---------------------------------
  # Preferred storage is G x J; also accept MixSim-style J x G and transpose.
  mu <- as.matrix(true_theta$mu)
  if (!is.numeric(mu) || anyNA(mu)) {
    stop("`true_theta$mu` must be numeric without missing values.")
  }
  if (nrow(mu) == G && ncol(mu) == n_celltypes) {
    mu_gj <- mu
  } else if (nrow(mu) == n_celltypes && ncol(mu) == G) {
    mu_gj <- t(mu)
  } else {
    stop("`true_theta$mu` must be G x J (or J x G).")
  }

  # --- Mixing proportions p: unit simplex ------------------------------------
  # Accept a single mixture weight vector (length J) or sample-wise matrices
  # (J x N or N x J). Matrix forms are reduced to length-J by averaging so that
  # MixSim / Jeffreys summaries still see one weight per cell type. Every
  # supplied simplex (vector, or each column/row of a matrix) must be
  # non-negative and sum to 1.
  p_vec <- NULL
  if (is.null(true_theta$p)) {
    if (isTRUE(require_p)) {
      stop("`true_theta` must contain `p`.", call. = FALSE)
    }
  } else {
    p_raw <- true_theta$p
    if (is.null(dim(p_raw))) {
      # Length-J vector on the simplex (shared mixture weights).
      p_vec <- as.numeric(p_raw)
      if (length(p_vec) != n_celltypes) {
        stop(
          "`true_theta$p` must have length J (= third dim of sigma/Theta).",
          call. = FALSE
        )
      }
    } else {
      p_mat <- as.matrix(p_raw)
      if (nrow(p_mat) == n_celltypes) {
        # J x N: each column is a sample-wise cellular ratio on the simplex;
        # average across samples for mixture-level summaries.
        if (
          any(abs(colSums(p_mat) - 1) > 100 * .Machine$double.eps) ||
            any(p_mat < 0) ||
            anyNA(p_mat)
        ) {
          stop(
            "`true_theta$p` as J x N must have non-negative columns ",
            "summing to 1.",
            call. = FALSE
          )
        }
        p_vec <- rowMeans(p_mat)
      } else if (ncol(p_mat) == n_celltypes) {
        # N x J: each row is a sample-wise ratio; average across samples.
        if (
          any(abs(rowSums(p_mat) - 1) > 100 * .Machine$double.eps) ||
            any(p_mat < 0) ||
            anyNA(p_mat)
        ) {
          stop(
            "`true_theta$p` as N x J must have non-negative rows ",
            "summing to 1.",
            call. = FALSE
          )
        }
        p_vec <- colMeans(p_mat)
      } else {
        stop(
          "`true_theta$p` must be length J, J x N, or N x J.",
          call. = FALSE
        )
      }
    }
    # Final simplex check on the (possibly averaged) length-J weight vector.
    if (
      any(p_vec < 0) ||
        anyNA(p_vec) ||
        abs(sum(p_vec) - 1) > 100 * .Machine$double.eps
    ) {
      stop(
        "`true_theta$p` must be non-negative and sum to 1 ",
        "(after averaging over samples if matrix).",
        call. = FALSE
      )
    }
  }

  # --- Cell-type count J -----------------------------------------------------
  # Infer from the third dim of sigma/Theta, or verify a caller-supplied J.
  if (is.null(J)) {
    J <- n_celltypes
  } else {
    J <- as.integer(J)
    if (length(J) != 1L || is.na(J) || J != n_celltypes) {
      stop(
        "`J` must match the third dimension of sigma/Theta.",
        call. = FALSE
      )
    }
  }
  if (J < 2L) {
    stop("`J` must be at least 2.", call. = FALSE)
  }

  # Standardised fields for downstream metrics (mu always G x J).
  list(
    p = p_vec,
    mu = mu_gj,
    sigma = if (has_sigma) true_theta$sigma else NULL,
    Theta = if (has_theta) true_theta$Theta else NULL,
    G = G,
    J = J
  )
}

#' Compute deconvolution benchmark metrics
#'
#' @description
#' Returns three blocks: composition / regression scores, Monte Carlo
#' summaries (ADEMP-style), and optimisation / runtime diagnostics. When
#' \eqn{\boldsymbol{p}^{\star}} is a matrix \eqn{J\times N} (or
#' \eqn{\hat{\boldsymbol{p}}} is), cell-type Pearson correlation and Monte
#' Carlo summaries are computed **across samples**, not within a single
#' composition.
#'
#' @inheritParams deconvolute_ratios_Marquardt_Levenberg
#' @param estimated_p Estimated proportions: a length-\eqn{J} vector or a
#'   numeric matrix \eqn{J\times N}.
#' @param true_ratios Optional ground truth with the same layout. When
#'   omitted, reconstitution scores compare
#'   \eqn{\hat{\boldsymbol{y}}=\boldsymbol{\mu}\hat{\boldsymbol{p}}} with
#'   \eqn{\boldsymbol{y}}.
#' @param se Optional standard errors matching `estimated_p`.
#' @param lower,upper Optional interval bounds matching `estimated_p`
#'   (used for coverage and mean width). When omitted but `se` is
#'   supplied, Wald intervals at `level` are used.
#' @param elapsed,memory_bytes,kkt_residual Optional per-sample runtime
#'   diagnostics (recycled to \eqn{N}).
#' @param numerical_converged Optional logical per sample: the solver
#'   returned a finite simplex vector.
#' @param loglik_hat,loglik_true Optional per-sample log-likelihoods used
#'   for theoretical convergence (regret
#'   \eqn{\ell(\boldsymbol{p}^{\star})-\ell(\hat{\boldsymbol{p}})}).
#' @param presence_threshold Threshold \eqn{\varepsilon} for presence /
#'   F1 / false-positive mass (default `1e-4`).
#' @param level Wald coverage level when `lower` / `upper` are omitted
#'   (default `0.95`).
#' @param coverage_interval Interval for the Monte Carlo coverage *rate*
#'   \eqn{\hat\pi}: `"wilson"` (default), `"wald"`, or `"agresti_coull"`.
#'   This is not the interval for \eqn{p_j} itself.
#' @param algorithm Character label stored on each metrics row (typically
#'   the solver name). Recycled to length one.
#'
#' @return A named list with:
#' * `regression`: `global` (one row per sample) and `cell_type` (one
#'   row per cell type);
#' * `monte_carlo`: one row per cell type;
#' * `optimisation`: one row per sample.
#'
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2,
#'              dimnames = list(paste0("g", 1:2), paste0("ct", 1:2)))
#' y <- drop(mu %*% c(0.4, 0.6))
#' compute_benchmark_metrics(y, mu, estimated_p = c(0.45, 0.55),
#'                           true_ratios = c(0.4, 0.6))
#' @importFrom rlang .data
#' @export
compute_benchmark_metrics <- function(
  y,
  mean_signature_matrix,
  estimated_p,
  true_ratios = NULL,
  se = NULL,
  lower = NULL,
  upper = NULL,
  elapsed = NULL,
  memory_bytes = NULL,
  kkt_residual = NULL,
  numerical_converged = NULL,
  loglik_hat = NULL,
  loglik_true = NULL,
  presence_threshold = 1e-4,
  level = 0.95,
  coverage_interval = "wilson",
  algorithm = NA_character_
) {
  cell_names <- colnames(mean_signature_matrix)
  if (is.null(cell_names)) {
    j_dim <- if (is.null(dim(estimated_p))) {
      length(estimated_p)
    } else {
      nrow(as.matrix(estimated_p))
    }
    cell_names <- paste0("ct", seq_len(j_dim))
  }
  p_hat <- .as_p_matrix(estimated_p, cell_names)
  n_samples <- ncol(p_hat)
  n_celltypes <- nrow(p_hat)
  p_true <- if (is.null(true_ratios)) {
    NULL
  } else {
    .as_p_matrix(true_ratios, cell_names, n_samples)
  }
  se_mat <- .as_optional_matrix(se, n_celltypes, n_samples, cell_names)
  bounds <- .wald_bounds(p_hat, se_mat, lower, upper, level)
  elapsed <- .recycle_numeric(elapsed, n_samples)
  memory_bytes <- .recycle_numeric(memory_bytes, n_samples)
  kkt_residual <- .recycle_numeric(kkt_residual, n_samples)
  loglik_hat <- .recycle_numeric(loglik_hat, n_samples)
  loglik_true <- .recycle_numeric(loglik_true, n_samples)
  if (is.null(numerical_converged)) {
    numerical_converged <- apply(p_hat, 2L, function(col) {
      all(is.finite(col))
    })
  }
  numerical_converged <- as.logical(.recycle_numeric(
    as.numeric(numerical_converged),
    n_samples
  ))
  regret <- loglik_true - loglik_hat
  theoretical_converged <- rep(NA, n_samples)
  finite_regret <- is.finite(regret)
  theoretical_converged[finite_regret] <- regret[finite_regret] <= 1e-3

  if (is.null(dim(y))) {
    y_mat <- matrix(as.numeric(y), ncol = 1L)
  } else {
    y_mat <- as.matrix(y)
  }
  if (ncol(y_mat) == 1L && n_samples > 1L) {
    y_mat <- y_mat[, rep(1L, n_samples), drop = FALSE]
  }

  global_rows <- lapply(seq_len(n_samples), function(i) {
    .regression_global_row(
      p_true = if (is.null(p_true)) NULL else p_true[, i],
      p_hat = p_hat[, i],
      y = y_mat[, i],
      mean_signature_matrix = mean_signature_matrix,
      sample_id = paste0("sample_", i),
      algorithm = algorithm
    )
  })
  regression_global <- dplyr::bind_rows(global_rows)

  cell_type <- .regression_cell_type_table(
    p_true = p_true,
    p_hat = p_hat,
    presence_threshold = presence_threshold,
    algorithm = algorithm
  )
  monte_carlo <- .monte_carlo_table(
    p_true = p_true,
    p_hat = p_hat,
    se = se_mat,
    lower = bounds$lower,
    upper = bounds$upper,
    coverage_interval = coverage_interval,
    algorithm = algorithm
  )

  opt <- tibble::tibble(
    sample_id = paste0("sample_", seq_len(n_samples)),
    algorithm = algorithm,
    elapsed_sec = elapsed,
    memory_bytes = memory_bytes,
    kkt_residual = kkt_residual,
    numerical_converged = numerical_converged,
    theoretical_converged = theoretical_converged,
    loglik_regret = regret
  )
  p_hat_tbl <- tibble::as_tibble(t(p_hat))
  optimisation <- dplyr::bind_cols(opt, p_hat_tbl)

  list(
    regression = list(global = regression_global, cell_type = cell_type),
    monte_carlo = monte_carlo,
    optimisation = optimisation
  )
}

#' @keywords internal
#' @noRd
.as_p_matrix <- function(p, cell_names, n_samples = NULL) {
  j <- length(cell_names)
  if (is.null(dim(p))) {
    p <- matrix(as.numeric(p), ncol = 1L)
  } else {
    p <- as.matrix(p)
    # Prefer J x N (cell types in rows). Transpose N x J only.
    if (nrow(p) != j && ncol(p) == j) {
      p <- t(p)
    }
  }
  if (nrow(p) != j) {
    stop(
      "`estimated_p` / `true_ratios` must have J rows (cell types).",
      call. = FALSE
    )
  }
  if (!is.null(n_samples) && ncol(p) == 1L && n_samples > 1L) {
    p <- p[, rep(1L, n_samples), drop = FALSE]
  }
  rownames(p) <- cell_names
  p
}

#' @keywords internal
#' @noRd
.as_optional_matrix <- function(x, n_celltypes, n_samples, cell_names) {
  if (is.null(x)) {
    out <- matrix(NA_real_, n_celltypes, n_samples)
    rownames(out) <- cell_names
    return(out)
  }
  .as_p_matrix(x, cell_names, n_samples)
}

#' @keywords internal
#' @noRd
.recycle_numeric <- function(x, n) {
  if (is.null(x)) {
    return(rep(NA_real_, n))
  }
  x <- as.numeric(x)
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop("Optional per-sample vectors must have length N.", call. = FALSE)
  }
  x
}

#' @keywords internal
#' @noRd
.wald_bounds <- function(p_hat, se, lower, upper, level) {
  n_samples <- ncol(p_hat)
  cell_names <- rownames(p_hat)
  if (!is.null(lower) && !is.null(upper)) {
    return(list(
      lower = .as_p_matrix(lower, cell_names, n_samples),
      upper = .as_p_matrix(upper, cell_names, n_samples)
    ))
  }
  z <- stats::qnorm((1 + level) / 2)
  list(
    lower = p_hat - z * se,
    upper = p_hat + z * se
  )
}

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

#' @keywords internal
#' @noRd
.regression_cell_type_table <- function(
  p_true,
  p_hat,
  presence_threshold,
  algorithm
) {
  cell_names <- rownames(p_hat)
  if (is.null(p_true)) {
    return(tibble::tibble(
      algorithm = character(),
      cell_type = character(),
      pearson = numeric(),
      presence_f1 = numeric(),
      false_positive_mass = numeric()
    ))
  }
  n_celltypes <- nrow(p_hat)
  n_samples <- ncol(p_hat)
  purrr::map_dfr(seq_len(n_celltypes), function(j) {
    counts <- list(tp = 0L, fp = 0L, fn = 0L, false_positive_mass = 0)
    for (i in seq_len(n_samples)) {
      one <- .presence_counts(
        p_true[j, i],
        p_hat[j, i],
        threshold = presence_threshold
      )
      counts$tp <- counts$tp + as.integer(one$tp)
      counts$fp <- counts$fp + as.integer(one$fp)
      counts$fn <- counts$fn + as.integer(one$fn)
      counts$false_positive_mass <- counts$false_positive_mass +
        one$false_positive_mass
    }
    tibble::tibble(
      algorithm = algorithm,
      cell_type = cell_names[[j]],
      pearson = .pearson_safe(p_true[j, ], p_hat[j, ]),
      presence_f1 = .f1_from_counts(counts$tp, counts$fp, counts$fn),
      false_positive_mass = counts$false_positive_mass / n_samples
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

#' Parallel deconvolution of a bulk expression matrix
#'
#' @description
#' For each column \eqn{\boldsymbol{y}_{\cdot i}} of the bulk matrix
#' \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}}, estimates
#' \eqn{\hat{\boldsymbol{p}}_{\cdot i}} with every supplied algorithm.
#' Samples are iterated with `furrr` when `cores > 1`; algorithms are
#' sequential, so workers are never nested. When covariance information
#' is provided, DeCovarT methods maximise
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})} under
#' \eqn{\boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})\sim
#' \mathcal{N}_{G}(
#'   \boldsymbol{\mu}\boldsymbol{p},
#'   \boldsymbol{\Sigma}(\boldsymbol{p})
#' )}.
#'
#' @inheritParams fit_decovart
#' @param true_ratios Optional ground-truth proportions: a length-\eqn{J}
#'   numeric vector (recycled over samples) or a numeric matrix
#'   \eqn{J\times N} (or \eqn{N\times J}).
#' @param Sigma Optional array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}}
#'   of cell-type covariances (numeric; off-diagonals may be negative).
#' @param deconvolution_functions Named list; each element has `FUN` and
#'   optional `additional_parameters`.
#' @param cores Number of `furrr` workers for the **sample** loop.
#'   Defaults to `getOption("mc.cores", 1L)`. Use `cores = 1` to stay
#'   sequential (CRAN examples and nested callers). With `cores > 1`,
#'   `.map_samples()` uses `furrr_options(seed = TRUE)`, which assigns
#'   an independent L'Ecuyer-CMRG stream per worker (R 4.1 / `future`
#'   streams). Do not offset a shared seed by the sample index.
#' @param coverage_interval Passed to [compute_benchmark_metrics()]
#'   (`"wilson"`, `"wald"`, or `"agresti_coull"`).
#' @param verbose If `TRUE`, emit progress for sequential sample loops
#'   (every `progress_every` samples). Parallel `furrr` runs log a
#'   one-line notice instead of per-sample ticks.
#' @param progress_every Sample-progress interval when `verbose` is
#'   `TRUE` and `cores = 1`. Defaults to `10L`.
#'
#' @return A named list from [compute_benchmark_metrics()] with
#'   `regression` (global and cell-type subtables), `monte_carlo`, and
#'   `optimisation` (per-sample elapsed time, memory, KKT residual, and
#'   \eqn{\hat{\boldsymbol{p}}}). First-generation solvers still call
#'   [repair_simplex()]; the three native DeCovarT maps (ILR or
#'   \eqn{p/\sum p}) already lie on the simplex.
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' Y <- simulate_bulk_mixture(mu, Sigma, p = c(0.5, 0.5), n = 2)$Y
#' deconvolute_ratios(
#'   signature_matrix = mu,
#'   bulk_expression = Y,
#'   true_ratios = c(0.5, 0.5),
#'   Sigma = Sigma,
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   cores = 1
#' )
#' @srrstats {G1.0} Closest published statistical method is DSection
#'   (Erkkila et al. 2010): univariate Gaussian convolution. DeMix, DeMixT
#'   and ISOpureR also treat cell-type profiles as latent; Ogundijo and Wang
#'   (2017) extend DSection by sequential Monte Carlo. DeCovarT generalises
#'   DSection to a multivariate convolution with sparse cell-type covariance.
#' @srrstats {G1.1} First implementation (in R or otherwise) of that
#'   multivariate Gaussian-convolution MLE for bulk deconvolution.
#' @srrstats {G2.8} Unique pre-processing gate:
#'   `.prepare_deconvolution_inputs()` then passes numeric matrices /
#'   arrays to every solver.
#' @srrstats {G2.13} Missing values in `y`, the signature, `Sigma` or
#'   `true_ratios` raise an error here, before any solver is called.
#' @srrstats {G2.14a} The missing-data policy is error-only (no imputation).
#' @srrstats {G2.15} After this check, `mean` / `cor` / `var` never receive
#'   incomplete expression data (default `na.rm = FALSE` is then safe).
#' @srrstats {G2.16} `NaN` / `Inf` / `-Inf` are rejected with the same error.
#' @srrstats {G5.8d} \eqn{J > G} raises an undetermined-mixture error.
#' @references
#' \insertRef{erkkilaProbabilisticAnalysisGene2010}{DeCovarT}
#' \insertRef{ahnDeMixDeconvolutionMixed2013}{DeCovarT}
#' \insertRef{wangTranscriptomeDeconvolutionHeterogeneous2018}{DeCovarT}
#' \insertRef{anghelISOpureRImplementationComputational2015}{DeCovarT}
#' \insertRef{ogundijoSequentialMonteCarlo2017}{DeCovarT}
#' \insertRef{hafemeisterNormalizationVarianceStabilization2019}{DeCovarT}
#' \insertRef{chionBayesianFrameworkMultivariate2023}{DeCovarT}
#' \insertRef{newmanRobustEnumerationCell2015}{DeCovarT}
#' @importFrom rlang .data
#' @export
#' @seealso [deconvolute_ratios_Marquardt_Levenberg()]
deconvolute_ratios <- function(
  signature_matrix,
  bulk_expression,
  true_ratios = NULL,
  Sigma = NULL,
  deconvolution_functions = NULL,
  standardise = FALSE,
  scaled = FALSE,
  cores = getOption("mc.cores", 1L),
  coverage_interval = "wilson",
  verbose = FALSE,
  progress_every = 10L
) {
  aligned <- .prepare_deconvolution_inputs(
    signature_matrix = signature_matrix,
    bulk_expression = bulk_expression,
    true_ratios = true_ratios,
    Sigma = Sigma,
    standardise = standardise,
    scaled = scaled
  )
  mean_signature_matrix <- aligned$signature_matrix
  Y <- aligned$bulk_expression
  true_ratios <- aligned$true_ratios
  Sigma <- aligned$Sigma
  cell_names <- colnames(mean_signature_matrix)
  n_celltypes <- ncol(mean_signature_matrix)

  if (
    is.null(deconvolution_functions) || length(deconvolution_functions) == 0L
  ) {
    .ui_abort("{.arg deconvolution_functions} must be a named list.")
  }

  per_algorithm <- purrr::imap(
    deconvolution_functions,
    function(deconvolution_function, algorithm) {
      additional_parameters <- deconvolution_function$additional_parameters
      sample_fits <- .map_samples(
        ncol(Y),
        function(i) {
          .fit_one_bulk_sample(
            i = i,
            Y = Y,
            true_ratios = true_ratios,
            mean_signature_matrix = mean_signature_matrix,
            Sigma = Sigma,
            deconvolution_function = deconvolution_function,
            additional_parameters = additional_parameters,
            cell_names = cell_names,
            n_celltypes = n_celltypes
          )
        },
        cores = cores,
        verbose = verbose,
        progress_every = progress_every,
        progress_label = algorithm
      )
      .assemble_algorithm_metrics(
        sample_fits = sample_fits,
        Y = Y,
        mean_signature_matrix = mean_signature_matrix,
        algorithm = algorithm,
        coverage_interval = coverage_interval
      )
    }
  )

  .bind_benchmark_lists(per_algorithm)
}

#' @keywords internal
#' @noRd
.fit_one_bulk_sample <- function(
  i,
  Y,
  true_ratios,
  mean_signature_matrix,
  Sigma,
  deconvolution_function,
  additional_parameters,
  cell_names,
  n_celltypes
) {
  y_i <- Y[, i, drop = TRUE]
  p_i <- if (is.null(true_ratios)) {
    NULL
  } else {
    true_ratios[, i, drop = TRUE]
  }
  list_arguments <- c(
    list(
      "y" = y_i,
      "mean_signature_matrix" = mean_signature_matrix,
      "Sigma" = Sigma
    ),
    additional_parameters
  )
  formal_args <- names(formals(deconvolution_function$FUN))
  mem0 <- .process_memory_bytes()
  t0 <- proc.time()[["elapsed"]]
  success_estimation <- tryCatch(
    {
      raw_out <- do.call(
        deconvolution_function$FUN,
        list_arguments[names(list_arguments) %in% formal_args]
      )
      estimated_p <- .coerce_estimated_p(raw_out, cell_names)
      se <- .extract_estimated_se(raw_out, n_celltypes)
      names(se) <- cell_names
      list(p = estimated_p, se = se, error = NULL)
    },
    error = function(e) {
      warning(conditionMessage(e), call. = FALSE)
      list(
        p = stats::setNames(rep(NA_real_, n_celltypes), cell_names),
        se = stats::setNames(rep(NA_real_, n_celltypes), cell_names),
        error = e
      )
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  mem1 <- .process_memory_bytes()
  memory_bytes <- max(c(mem0, mem1), na.rm = TRUE)
  if (!is.finite(memory_bytes)) {
    memory_bytes <- NA_real_
  }

  kkt <- NA_real_
  loglik_hat <- NA_real_
  loglik_true <- NA_real_
  p_hat <- success_estimation$p
  numerical_converged <- is.null(success_estimation$error) &&
    all(is.finite(p_hat))
  if (numerical_converged && !is.null(Sigma)) {
    kkt <- tryCatch(
      {
        grad <- gradient_loglik_unconstrained(
          p_hat,
          y_i,
          mean_signature_matrix,
          Sigma
        )
        .kkt_residual(p_hat, grad)
      },
      error = function(e) NA_real_
    )
    loglik_hat <- tryCatch(
      loglik_multivariate(p_hat, y_i, mean_signature_matrix, Sigma),
      error = function(e) NA_real_
    )
    if (!is.null(p_i)) {
      loglik_true <- tryCatch(
        loglik_multivariate(p_i, y_i, mean_signature_matrix, Sigma),
        error = function(e) NA_real_
      )
    }
  }

  list(
    p = p_hat,
    se = success_estimation$se,
    true_p = p_i,
    y = y_i,
    elapsed = elapsed,
    memory_bytes = memory_bytes,
    kkt_residual = kkt,
    numerical_converged = numerical_converged,
    loglik_hat = loglik_hat,
    loglik_true = loglik_true
  )
}

#' @keywords internal
#' @noRd
.assemble_algorithm_metrics <- function(
  sample_fits,
  Y,
  mean_signature_matrix,
  algorithm,
  coverage_interval = "wilson"
) {
  p_hat <- do.call(cbind, purrr::map(sample_fits, "p"))
  se <- do.call(cbind, purrr::map(sample_fits, "se"))
  true_list <- purrr::map(sample_fits, "true_p")
  true_ratios <- if (all(purrr::map_lgl(true_list, is.null))) {
    NULL
  } else {
    do.call(cbind, true_list)
  }
  compute_benchmark_metrics(
    y = Y,
    mean_signature_matrix = mean_signature_matrix,
    estimated_p = p_hat,
    true_ratios = true_ratios,
    se = se,
    elapsed = purrr::map_dbl(sample_fits, "elapsed"),
    memory_bytes = purrr::map_dbl(sample_fits, "memory_bytes"),
    kkt_residual = purrr::map_dbl(sample_fits, "kkt_residual"),
    numerical_converged = purrr::map_lgl(
      sample_fits,
      "numerical_converged"
    ),
    loglik_hat = purrr::map_dbl(sample_fits, "loglik_hat"),
    loglik_true = purrr::map_dbl(sample_fits, "loglik_true"),
    coverage_interval = coverage_interval,
    algorithm = algorithm
  )
}

#' @keywords internal
#' @noRd
.bind_benchmark_lists <- function(per_algorithm) {
  list(
    regression = list(
      global = dplyr::bind_rows(
        purrr::map(per_algorithm, \(x) x$regression$global)
      ),
      cell_type = dplyr::bind_rows(
        purrr::map(per_algorithm, \(x) x$regression$cell_type)
      )
    ),
    monte_carlo = dplyr::bind_rows(
      purrr::map(per_algorithm, "monte_carlo")
    ),
    optimisation = dplyr::bind_rows(
      purrr::map(per_algorithm, "optimisation")
    )
  )
}
