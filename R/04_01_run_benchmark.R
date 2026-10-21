# Benchmark driver: fits every algorithm to every bulk sample, records
# runtime and optimisation diagnostics, and assembles the global
# (04_03_performance_metrics_global.R) and cell-type
# (04_04_performance_metrics_celltype.R) performance tables.

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
#'   Convolution-likelihood solvers fill these from [vcov_ilr_delta()]
#'   (expected Fisher information pulled back through the ILR chart).
#'   Mean-only solvers (NNLS, LSEI) leave them missing.
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
#' @param converged Optional list of length \eqn{N} with homogenised
#'   optimiser fields (`iterations`, `criterion`, `code`, `rdm`).
#'   Closed-form solvers store `NA` entries.
#' @param local_maximum Optional logical per sample: ILR score small and
#'   Hessian negative definite (\eqn{\lambda_{\min}<0} and
#'   \eqn{\lambda_{\max}<0}).
#' @param min_eigenvalue,max_eigenvalue Optional ILR Hessian extrema.
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
#' * `monte_carlo`: one row per cell type, including `theoretical_se`
#'   (filled later from expected Fisher at \eqn{p^{\star}});
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
  converged = NULL,
  local_maximum = NULL,
  min_eigenvalue = NULL,
  max_eigenvalue = NULL,
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
  se_mat <- if (is.null(se)) {
    matrix(
      NA_real_,
      n_celltypes,
      n_samples,
      dimnames = list(cell_names, NULL)
    )
  } else {
    .as_p_matrix(se, cell_names, n_samples)
  }
  # Supplied bounds are used as they are; otherwise Wald intervals at `level`.
  if (!is.null(lower) && !is.null(upper)) {
    lower <- .as_p_matrix(lower, cell_names, n_samples)
    upper <- .as_p_matrix(upper, cell_names, n_samples)
  } else {
    z <- stats::qnorm((1 + level) / 2)
    lower <- p_hat - z * se_mat
    upper <- p_hat + z * se_mat
  }
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
  local_maximum <- if (is.null(local_maximum)) {
    rep(NA, n_samples)
  } else {
    as.logical(.recycle_numeric(as.numeric(local_maximum), n_samples))
  }
  min_eigenvalue <- .recycle_numeric(min_eigenvalue, n_samples)
  max_eigenvalue <- .recycle_numeric(max_eigenvalue, n_samples)
  if (is.null(converged)) {
    converged <- replicate(
      n_samples,
      .empty_convergence(),
      simplify = FALSE
    )
  }
  if (!is.list(converged) || length(converged) != n_samples) {
    stop("`converged` must be a list of length N.", call. = FALSE)
  }

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
    lower = lower,
    upper = upper,
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
    loglik_regret = regret,
    local_maximum = local_maximum,
    min_eigenvalue = min_eigenvalue,
    max_eigenvalue = max_eigenvalue,
    converged = converged
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
#'   \eqn{p/\sum p}) already lie on the simplex. When `Sigma` is
#'   supplied, Newton-Raphson, Marquardt-Levenberg, L-BFGS-B, BFGS
#'   (gradient), and simulated annealing attach ILR Wald standard
#'   errors from [vcov_ilr_delta()] so `monte_carlo` can report
#'   coverage, mean model SE, and the SE/SD ratio. Those SEs are the
#'   expected-Fisher delta-method bounds used by
#'   [confint.decovart_fit()], not diagonal entries of the observed
#'   Hessian. Mean-only solvers (NNLS, LSEI) leave Wald columns
#'   missing: those estimators do not use the convolution likelihood.
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

  if (
    is.null(deconvolution_functions) || length(deconvolution_functions) == 0L
  ) {
    .ui_abort("{.arg deconvolution_functions} must be a named list.")
  }

  per_algorithm <- purrr::imap(
    deconvolution_functions,
    function(deconvolution_function, algorithm) {
      fits <- .map_samples(
        ncol(Y),
        function(i) {
          .fit_one_bulk_sample(
            y_i = Y[, i, drop = TRUE],
            p_i = if (is.null(true_ratios)) NULL else true_ratios[, i],
            mean_signature_matrix = mean_signature_matrix,
            Sigma = Sigma,
            deconvolution_function = deconvolution_function,
            cell_names = cell_names
          )
        },
        cores = cores,
        verbose = verbose,
        progress_every = progress_every,
        progress_label = algorithm
      )
      true_list <- purrr::map(fits, "true_p")
      compute_benchmark_metrics(
        y = Y,
        mean_signature_matrix = mean_signature_matrix,
        estimated_p = do.call(cbind, purrr::map(fits, "p")),
        true_ratios = if (all(purrr::map_lgl(true_list, is.null))) {
          NULL
        } else {
          do.call(cbind, true_list)
        },
        se = do.call(cbind, purrr::map(fits, "se")),
        elapsed = purrr::map_dbl(fits, "elapsed"),
        memory_bytes = purrr::map_dbl(fits, "memory_bytes"),
        kkt_residual = purrr::map_dbl(fits, "kkt_residual"),
        numerical_converged = purrr::map_lgl(fits, "numerical_converged"),
        loglik_hat = purrr::map_dbl(fits, "loglik_hat"),
        loglik_true = purrr::map_dbl(fits, "loglik_true"),
        converged = purrr::map(fits, "converged"),
        local_maximum = purrr::map_lgl(fits, "local_maximum"),
        min_eigenvalue = purrr::map_dbl(fits, "min_eigenvalue"),
        max_eigenvalue = purrr::map_dbl(fits, "max_eigenvalue"),
        coverage_interval = coverage_interval,
        algorithm = algorithm
      )
    }
  )
  bind <- function(get) dplyr::bind_rows(purrr::map(per_algorithm, get))
  list(
    regression = list(
      global = bind(\(x) x$regression$global),
      cell_type = bind(\(x) x$regression$cell_type)
    ),
    monte_carlo = bind("monte_carlo"),
    optimisation = bind("optimisation")
  )
}

#' Fit one bulk sample with one algorithm and record its diagnostics
#'
#' Convolution-likelihood solvers get ILR Wald standard errors from the
#' expected Fisher information; mean-only solvers leave them missing. When
#' `Sigma` is supplied, the KKT residual, the log-likelihoods at
#' \eqn{\hat{\boldsymbol{p}}} and \eqn{\boldsymbol{p}^{\star}}, and the ILR
#' Hessian extrema are recorded as well.
#'
#' @keywords internal
#' @noRd
.fit_one_bulk_sample <- function(
  y_i,
  p_i,
  mean_signature_matrix,
  Sigma,
  deconvolution_function,
  cell_names
) {
  n_celltypes <- length(cell_names)
  fun <- deconvolution_function$FUN
  list_arguments <- c(
    list(y = y_i, mean_signature_matrix = mean_signature_matrix, Sigma = Sigma),
    deconvolution_function$additional_parameters
  )
  formal_args <- names(formals(fun))
  if ("return_model" %in% formal_args) {
    list_arguments$return_model <- TRUE
  }
  convolution_solver <- any(vapply(
    list(
      deconvolute_ratios_Newton_Raphson,
      deconvolute_ratios_Marquardt_Levenberg,
      deconvolute_ratios_L_BFGS_B,
      deconvolute_ratios_gradient_descent,
      deconvolute_ratios_simulated_annealing
    ),
    identical,
    logical(1),
    fun
  ))
  na_vector <- stats::setNames(rep(NA_real_, n_celltypes), cell_names)
  mem0 <- .process_memory_bytes()
  t0 <- proc.time()[["elapsed"]]
  fit <- tryCatch(
    {
      raw_out <- do.call(
        fun,
        list_arguments[names(list_arguments) %in% formal_args]
      )
      estimated <- if (is.list(raw_out) && !is.null(raw_out$coefficients)) {
        raw_out$coefficients
      } else {
        raw_out
      }
      estimated_p <- stats::setNames(as.numeric(estimated), cell_names)
      se <- if (
        convolution_solver && !is.null(Sigma) && all(is.finite(estimated_p))
      ) {
        .ilr_wald_se(estimated_p, mean_signature_matrix, Sigma, warn = FALSE)
      } else {
        na_vector
      }
      list(
        p = estimated_p,
        se = se,
        ok = TRUE,
        converged = .homogenise_convergence(raw_out)
      )
    },
    error = function(e) {
      warning(conditionMessage(e), call. = FALSE)
      list(
        p = na_vector,
        se = na_vector,
        ok = FALSE,
        converged = .empty_convergence()
      )
    }
  )
  elapsed <- proc.time()[["elapsed"]] - t0
  memory_bytes <- suppressWarnings(max(
    c(mem0, .process_memory_bytes()),
    na.rm = TRUE
  ))
  if (!is.finite(memory_bytes)) {
    memory_bytes <- NA_real_
  }

  out <- list(
    p = fit$p,
    se = fit$se,
    true_p = p_i,
    elapsed = elapsed,
    memory_bytes = memory_bytes,
    kkt_residual = NA_real_,
    numerical_converged = fit$ok && all(is.finite(fit$p)),
    loglik_hat = NA_real_,
    loglik_true = NA_real_,
    converged = fit$converged,
    local_maximum = NA,
    min_eigenvalue = NA_real_,
    max_eigenvalue = NA_real_
  )
  if (out$numerical_converged && !is.null(Sigma)) {
    safely <- function(expr) tryCatch(expr, error = function(e) NA_real_)
    p_hat <- fit$p
    out$kkt_residual <- safely(.kkt_residual(
      p_hat,
      gradient_loglik_unconstrained(p_hat, y_i, mean_signature_matrix, Sigma)
    ))
    out$loglik_hat <- safely(
      loglik_multivariate(p_hat, y_i, mean_signature_matrix, Sigma)
    )
    if (!is.null(p_i)) {
      out$loglik_true <- safely(
        loglik_multivariate(p_i, y_i, mean_signature_matrix, Sigma)
      )
    }
    diagnostics <- tryCatch(
      boundary_diagnostics(p_hat, y_i, mean_signature_matrix, Sigma),
      error = function(e) NULL
    )
    if (!is.null(diagnostics)) {
      out$local_maximum <- isTRUE(diagnostics$local_maximum)
      out$min_eigenvalue <- diagnostics$min_eigenvalue
      out$max_eigenvalue <- diagnostics$max_eigenvalue
    }
  }
  out
}

#' Map over bulk samples, optionally with `furrr`
#'
#' @keywords internal
#' @noRd
.map_samples <- function(
  n,
  fun,
  cores,
  verbose = FALSE,
  progress_every = 10L,
  progress_label = NULL
) {
  idx <- seq_len(n)
  cores <- as.integer(cores)
  if (length(cores) != 1L || is.na(cores) || cores < 1L) {
    .ui_abort("{.arg cores} must be a positive integer.")
  }
  if (cores == 1L) {
    if (!isTRUE(verbose) || n < 1L) {
      return(lapply(idx, fun))
    }
    every <- as.integer(progress_every)
    if (length(every) != 1L || is.na(every) || every < 1L) {
      every <- 10L
    }
    prefix <- if (is.null(progress_label) || !nzchar(progress_label)) {
      ""
    } else {
      paste0(progress_label, " \u00b7 ")
    }
    return(lapply(idx, function(i) {
      out <- fun(i)
      if (i == n || i == 1L || (i %% every == 0L)) {
        .ui_info("{prefix}inferred sample {.val {i}}/{.val {n}}.")
      }
      out
    }))
  }
  if (isTRUE(verbose)) {
    .ui_info(
      "Mapping {.val {n}} samples on {.val {cores}} workers (no per-sample ticks)."
    )
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

#' Process memory in bytes (PSS on Linux when `ps` is available)
#'
#' Windows `ps` objects have no `pss` field; `x[["pss"]]` on a named
#' numeric vector then errors with subscript-out-of-bounds.
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
  .one_named <- function(nm) {
    if (!nm %in% names(info)) {
      return(NULL)
    }
    val <- info[[nm]]
    if (!is.null(val) && length(val) == 1L && is.finite(val)) {
      return(as.numeric(val))
    }
    NULL
  }
  pss <- .one_named("pss")
  if (!is.null(pss)) {
    return(pss)
  }
  rss <- .one_named("rss")
  if (!is.null(rss)) {
    return(rss)
  }
  NA_real_
}

#' Empty homogenised optimiser record
#'
#' @keywords internal
#' @noRd
.empty_convergence <- function() {
  list(
    iterations = NA_integer_,
    criterion = NA_character_,
    code = NA_integer_,
    rdm = NA_real_
  )
}

#' First finite scalar of the requested type
#'
#' @keywords internal
#' @noRd
.first_scalar <- function(..., type = "integer") {
  vals <- list(...)
  for (v in vals) {
    if (is.null(v) || length(v) < 1L) {
      next
    }
    x <- v[[1L]]
    if (identical(type, "character")) {
      x <- as.character(x)
      if (is.na(x) || !nzchar(x)) {
        next
      }
      return(x)
    }
    if (identical(type, "double")) {
      x <- suppressWarnings(as.numeric(x))
      if (length(x) == 1L && is.finite(x)) {
        return(x)
      }
      next
    }
    x <- suppressWarnings(as.integer(x))
    if (length(x) == 1L && !is.na(x)) {
      return(x)
    }
  }
  if (identical(type, "character")) {
    return(NA_character_)
  }
  if (identical(type, "double")) {
    return(NA_real_)
  }
  NA_integer_
}

#' Homogenised optimiser diagnostics (iterative solvers only)
#'
#' Closed-form methods (NNLS, LSEI, CIBERSORT) return a plain vector without
#' a `convergence` record, so their fields are `NA`.
#'
#' @keywords internal
#' @noRd
.homogenise_convergence <- function(raw_out) {
  conv <- if (is.list(raw_out)) raw_out$convergence else NULL
  if (!is.list(conv)) {
    return(.empty_convergence())
  }
  list(
    iterations = .first_scalar(conv$iterations, conv$ni, type = "integer"),
    criterion = .first_scalar(
      conv$criterion,
      conv$message,
      type = "character"
    ),
    code = .first_scalar(
      conv$code,
      conv$istop,
      conv$convergence,
      type = "integer"
    ),
    rdm = .first_scalar(conv$rdm, type = "double")
  )
}
