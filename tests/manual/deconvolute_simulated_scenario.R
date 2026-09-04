# ---- tests/manual/deconvolute_simulated_scenario.R --------------------------
##
## TEMPORARY solver smoke test (Rbuildignored). Not part of testthat.
## Simulate bulk mixtures from the variance-driven GRN produced by
## scripts/fig03_variance_driven.R, with imbalanced p = (0.4, 0.4, 0.2).
## Cell types 1 and 2 are near mean-collinear; recovering both masses
## rather than lumping them is a check of the covariance-aware signal.
##
## Prerequisite: run fig03 at least once so
## data/synthetic_networks/true_grn_moments.rds exists.
##
## Five solvers from R/03_03_DeCovarT_estimate_ratios_frequentist.R:
##   * deconvolute_ratios_simulated_annealing() -- stochastic search on rho
##     ([stats::optim()], `method = "SANN"`); no gradient used.
##   * deconvolute_ratios_L_BFGS_B()            -- box-constrained ascent
##     directly on p ([stats::optim()], `method = "L-BFGS-B"`); analytic
##     gradient, safeguarded against singular Sigma(p).
##   * deconvolute_ratios_gradient_descent()    -- BFGS quasi-Newton ascent
##     on rho using the ANALYTIC gradient gradient_loglik_constrained().
##   * deconvolute_ratios_Marquardt_Levenberg() -- second-order ascent on rho
##     ([marqLevAlg::marqLevAlg()]) using the ANALYTIC gradient AND Hessian
##     (correctly negated for marqLevAlg's minimisation-convention Cholesky
##     routines; see the in-code comment) and structured istop/rdm checks
##     per Commenges et al. (2006) rather than a second run + regex scrape.
##   * deconvolute_ratios_Newton_Raphson()      -- second-order ascent on rho
##     ([stats::nlminb()]) using the ANALYTIC gradient AND Hessian.
## Every (solver, bulk sample) pair is deconvolved independently and
## defensively: both R errors (e.g. `solve()` on a near-singular Sigma(p),
## `repair_simplex()` rejecting an estimate) and warnings (e.g. boundary
## underflow of p) are caught without aborting the run, a stalled optimiser
## (unmoved initial guess) is flagged explicitly, and wall-clock time per
## call is recorded to compare computational cost.

devtools::load_all(".", quiet = TRUE)
set.seed(20260807)

true_grn_moments <- readRDS(
  file.path("data", "synthetic_networks", "true_grn_moments.rds")
)

## Restrict to an NSGA-II panel when present (legacy cache); otherwise
## use the full G = 50 design from fig03.
final_panel <- true_grn_moments$final_panel
if (is.null(final_panel)) {
  final_panel <- rownames(true_grn_moments$mu)
}
mu <- true_grn_moments$mu[final_panel, , drop = FALSE]
Sigma <- true_grn_moments$Sigma[final_panel, final_panel, , drop = FALSE]
celltype_names <- colnames(mu)
n_celltypes <- length(celltype_names)

n_bulk <- 20L # number of bulk mixtures to simulate
## Imbalanced true composition: 0.4 / 0.4 for the two mean-collinear cell
## types (1 and 2), 0.2 for cell type 3 (its own marker genes).
p_true <- stats::setNames(c(0.4, 0.4, 0.2), celltype_names)

bulk_mixtures <- simulate_bulk_mixture(
  signature_matrix = mu,
  Sigma = Sigma,
  p = p_true,
  n = n_bulk
)

message(
  "Simulated ",
  n_bulk,
  " bulk mixtures on the ",
  length(final_panel),
  "-gene curated panel (p_true = ",
  paste(paste0(celltype_names, "=", signif(p_true, 3)), collapse = " / "),
  ")."
)

## ---- Gradient sanity check: analytic vs numDeriv --------------------------
## deconvolute_ratios_gradient_descent() relies on the analytic
## gradient_loglik_constrained(). Cross-check it against numDeriv::grad()
## on loglik_multivariate_constrained() at the TRUE rho for one bulk
## sample, mirroring tests/testthat/test-03_03_DeCovarT.R, before trusting
## the BFGS run below.
true_rho <- additive_log_ratio(p_true)
gradient_numerical <- numDeriv::grad(
  loglik_multivariate_constrained,
  true_rho,
  method = "Richardson",
  method.args = list(eps = 1e-4, r = 6),
  y = bulk_mixtures$Y[, 1L],
  mean_signature_matrix = mu,
  Sigma = Sigma
)
gradient_analytic <- as.numeric(gradient_loglik_constrained(
  true_rho,
  bulk_mixtures$Y[, 1L],
  mu,
  Sigma
))
gradient_max_abs_diff <- max(abs(gradient_numerical - gradient_analytic))
message(
  "Gradient sanity check (sample_1, true rho): analytic = ",
  paste(signif(gradient_analytic, 4), collapse = ", "),
  "; numDeriv = ",
  paste(signif(gradient_numerical, 4), collapse = ", "),
  "; max |diff| = ",
  signif(gradient_max_abs_diff, 3),
  if (gradient_max_abs_diff < 1e-4) {
    " (matches, < 1e-4)."
  } else {
    " (WARNING: exceeds 1e-4 tolerance)."
  }
)

## ---- Deconvolution, one (solver, bulk sample) pair at a time -------------
## Each solver only accepts a single bulk vector `y`; loop explicitly over
## the n_bulk columns of Y for every solver. All solvers start from the
## same equi-balanced guess `initial_p`; a solver returning that guess
## unmoved (0 accepted / effective moves within `itmax`) is flagged as
## "stalled" rather than mistaken for a genuine (and here, since p_true is
## imbalanced, clearly wrong) estimate.
initial_p <- rep(1 / n_celltypes, n_celltypes)

solvers <- list(
  simulated_annealing = deconvolute_ratios_simulated_annealing,
  L_BFGS_B = deconvolute_ratios_L_BFGS_B,
  gradient_descent = deconvolute_ratios_gradient_descent,
  Marquardt_Levenberg = deconvolute_ratios_Marquardt_Levenberg,
  Newton_Raphson = deconvolute_ratios_Newton_Raphson
)

## Wall-clock time (ms) is recorded around the solver call only (excludes
## tibble bookkeeping) so that per-algorithm cost is directly comparable.
deconvolve_one_sample <- function(solver_fn, y, mean_signature_matrix, Sigma) {
  warnings_seen <- character()
  elapsed_ms <- NA_real_
  estimated_p <- tryCatch(
    withCallingHandlers(
      {
        start_time <- proc.time()[["elapsed"]]
        out <- solver_fn(
          y = y,
          mean_signature_matrix = mean_signature_matrix,
          Sigma = Sigma
        )
        # Plain `<-`, not `<<-`: this block is evaluated directly in
        # deconvolve_one_sample()'s own frame (a `{ }` block passed as an
        # argument introduces no new scope), so `elapsed_ms` here already
        # *is* the local variable declared above -- `<<-` would instead
        # skip that frame and leak a stray binding into the global env.
        elapsed_ms <- 1000 * (proc.time()[["elapsed"]] - start_time)
        out
      },
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  list(
    estimate = estimated_p,
    warnings = warnings_seen,
    elapsed_ms = elapsed_ms
  )
}

## ---- Assemble one row per (algorithm, bulk sample), flagging issues ------
deconvolution_results <- purrr::imap_dfr(solvers, function(solver_fn, algo) {
  purrr::imap_dfr(seq_len(n_bulk), function(i, .idx) {
    run <- deconvolve_one_sample(
      solver_fn,
      y = bulk_mixtures$Y[, i],
      mean_signature_matrix = mu,
      Sigma = Sigma
    )
    estimated_p <- run$estimate
    n_warnings <- length(run$warnings)

    if (inherits(estimated_p, "error")) {
      row <- tibble::tibble(
        algorithm = algo,
        sample = paste0("sample_", i),
        status = "error",
        n_warnings = n_warnings,
        elapsed_ms = run$elapsed_ms,
        sum_p = NA_real_,
        max_abs_error = NA_real_,
        note = conditionMessage(estimated_p)
      )
      for (ct in celltype_names) {
        row[[ct]] <- NA_real_
      }
      return(row)
    }

    sum_p <- sum(estimated_p)
    issues <- character()
    if (abs(sum_p - 1) > 1e-6) {
      issues <- c(issues, sprintf("sum(p) = %.6f", sum_p))
    }
    if (any(estimated_p < 0 | estimated_p > 1)) {
      issues <- c(issues, "estimate outside [0, 1]")
    }
    if (n_warnings > 0L) {
      issues <- c(issues, sprintf("%d optim warning(s)", n_warnings))
    }
    if (
      isTRUE(all.equal(
        as.numeric(estimated_p[celltype_names]),
        initial_p,
        tolerance = 1e-8
      ))
    ) {
      issues <- c(
        issues,
        "stalled at initial guess (0 accepted/effective moves)"
      )
    }

    tibble::tibble(
      algorithm = algo,
      sample = paste0("sample_", i),
      status = if (length(issues)) "flagged" else "ok",
      n_warnings = n_warnings,
      elapsed_ms = run$elapsed_ms,
      sum_p = sum_p,
      max_abs_error = max(abs(estimated_p[celltype_names] - p_true)),
      note = paste(issues, collapse = "; ")
    ) |>
      dplyr::bind_cols(tibble::as_tibble_row(estimated_p[celltype_names]))
  })
})

## ---- Per-algorithm summary -------------------------------------------------
algorithm_summary <- deconvolution_results |>
  dplyr::group_by(.data$algorithm) |>
  dplyr::summarise(
    n_ok = sum(.data$status == "ok"),
    n_flagged = sum(.data$status == "flagged"),
    n_errors = sum(.data$status == "error"),
    mean_max_abs_error = mean(.data$max_abs_error, na.rm = TRUE),
    mean_elapsed_ms = mean(.data$elapsed_ms, na.rm = TRUE),
    .groups = "drop"
  )

for (i in seq_len(nrow(algorithm_summary))) {
  row <- algorithm_summary[i, ]
  message(
    row$algorithm,
    " over ",
    n_bulk,
    " bulk samples: ",
    row$n_ok,
    " ok, ",
    row$n_flagged,
    " flagged, ",
    row$n_errors,
    " failed; mean max |error| = ",
    signif(row$mean_max_abs_error, 3),
    "; mean time = ",
    signif(row$mean_elapsed_ms, 3),
    " ms/sample."
  )
}
n_flagged_total <- sum(deconvolution_results$status != "ok")
if (n_flagged_total > 0L) {
  message(
    "Numerical inconsistencies detected across solvers -- see the `note` ",
    "column (e.g. stalled search, boundary p underflow, sum(p) drift, or ",
    "an istop/rdm warning from marqLevAlg::marqLevAlg())."
  )
}

## ---- Report ---------------------------------------------------------------
summary_caption <- paste0(
  "Per-algorithm summary over ",
  n_bulk,
  " synthetic bulk mixtures (true p = ",
  paste(paste0(celltype_names, "=", p_true), collapse = ", "),
  "), on the NSGA-II-curated ",
  length(final_panel),
  "-gene panel."
)
print(
  tinytable::tt(algorithm_summary, digits = 3, caption = summary_caption) |>
    tinytable::style_tt(j = 2:6, align = "c")
)

detail_caption <- paste0(
  "Per-sample deconvolution detail (",
  paste(names(solvers), collapse = ", "),
  "). ",
  if (n_flagged_total > 0L) {
    paste0(
      n_flagged_total,
      " (algorithm, sample) pair(s) were flagged or failed; see `note`."
    )
  } else {
    "No numerical inconsistencies were detected."
  }
)
numeric_cols <- c(
  "n_warnings",
  "elapsed_ms",
  "sum_p",
  "max_abs_error",
  celltype_names
)
print(
  tinytable::tt(deconvolution_results, digits = 3, caption = detail_caption) |>
    tinytable::style_tt(
      j = match(numeric_cols, names(deconvolution_results)),
      align = "c"
    )
)
