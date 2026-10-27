#!/usr/bin/env Rscript
# Numerical checks of the confidence-interval menus of both models.
#
#   1. Information matrices against independent computations:
#      population expected Fisher information against the Monte Carlo mean
#      of the negative observed Hessian; single-cell profile and nuisance
#      covariances against numDeriv Hessians of the full log-likelihood
#      in (rho, c, phi).
#   2. Monte Carlo coverage of the interior intervals:
#      population convolution (N pooled columns, J = 2 and J = 3) with
#      methods ilr, wald, firth, profile; single-cell mixture (the four
#      scenarios of single_cell_scenario_design(), a = 1) with methods
#      ilr, wald, profile, plus the intervals for c and phi.
#   3. Bootstrap spot checks (parametric and reference, fewer replicates).
#   4. Size of the boundary likelihood-ratio tests at a face or a vertex.
#
# Usage (from the repository root):
#   Rscript scripts/benchmark_metrics/check_confidence_interval_coverage.R [n_rep]
#
# Output: output/ci_checks/ci_checks.rds and a printed summary.

devtools::load_all(".", quiet = TRUE)

args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 1L) as.integer(args[[1L]]) else 300L
n_rep_boot <- max(20L, n_rep %/% 6L)
level <- 0.95
seed <- 20261022L
output_dir <- file.path("output", "ci_checks")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

mc_se <- function(rate, n) sqrt(rate * (1 - rate) / n)
covers <- function(interval, truth) {
  interval[, 1L] <= truth & truth <= interval[, 2L]
}
timed <- function(label, expr) {
  t0 <- proc.time()[["elapsed"]]
  out <- expr
  message(sprintf("%-40s %6.1f s", label, proc.time()[["elapsed"]] - t0))
  out
}

# --- 1. Information matrices ---------------------------------------------

population_moments <- function(n_celltypes, n_genes = 12L) {
  withr::with_seed(seed + n_celltypes, {
    simulate_hierarchical_grn_moments(
      n_genes = n_genes,
      n_celltypes = n_celltypes,
      target_cosine = 0.3,
      precision_shift = 0.3,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
  })
}

check_population_information <- function(n_celltypes, p, n_draws = 4000L) {
  moments <- population_moments(n_celltypes)
  mu <- moments$mean_profiles
  sigma <- moments$covariance_matrices
  y <- withr::with_seed(seed, {
    simulate_bulk_mixture(mu, sigma, p = p, n = n_draws)$Y
  })
  observed <- Reduce(
    `+`,
    lapply(seq_len(n_draws), function(i) {
      -hessian_loglik_unconstrained(p, y[, i], mu, sigma)
    })
  ) /
    n_draws
  expected <- expected_fisher_unconstrained(p, mu, sigma)
  data.frame(
    check = sprintf("population I(p), J = %d", n_celltypes),
    max_relative_error = max(abs(observed - expected)) / max(abs(expected)),
    reference = sprintf("Monte Carlo mean of -H over %d draws", n_draws)
  )
}

single_cell_reference <- function(n_celltypes, cell_size) {
  moments <- withr::with_seed(seed + 10L * n_celltypes, {
    simulate_hierarchical_grn_moments(
      n_genes = 10L,
      n_celltypes = n_celltypes,
      mean_scale = 10,
      target_cosine = 0,
      nonnegative = TRUE,
      precision_shift = 0.5,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
  })
  rescale_single_cell_moments(moments, cell_size)
}

check_single_cell_information <- function(n_celltypes) {
  ref <- single_cell_reference(n_celltypes, seq_len(n_celltypes))
  p <- seq(n_celltypes, 1) / sum(seq_len(n_celltypes))
  sim <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = p,
    n_cells = 1e4,
    platform_scale = 1,
    n_samples = 1L,
    seed = seed
  )
  y <- sim$Y[, 1L]
  mu <- ref$mean_signature_matrix
  sigma <- ref$Sigma
  fit <- deconvolute_single_cell(y, mu, sigma)
  theta <- c(fit$rho, log(fit$mean_scale), log(fit$variance_scale))
  n_rho <- length(fit$rho)
  full_loglik <- function(theta) {
    loglik_single_cell(
      isometric_logistic(theta[seq_len(n_rho)]),
      y,
      mu,
      sigma,
      exp(theta[[n_rho + 1L]]),
      exp(theta[[n_rho + 2L]])
    )
  }
  joint <- solve(-numDeriv::hessian(full_loglik, theta))
  # Delta method from (log c, log phi) to (c, phi).
  scale_jac <- diag(c(fit$mean_scale, fit$variance_scale))
  nuisance <- scale_jac %*%
    joint[n_rho + 1:2, n_rho + 1:2, drop = FALSE] %*%
    scale_jac
  rel <- function(a, b) max(abs(a - b)) / max(abs(b))
  data.frame(
    check = sprintf(
      c("single-cell vcov_rho, J = %d", "single-cell vcov_nuisance, J = %d"),
      n_celltypes
    ),
    max_relative_error = c(
      rel(fit$vcov_rho, joint[seq_len(n_rho), seq_len(n_rho), drop = FALSE]),
      rel(fit$vcov_nuisance, nuisance)
    ),
    reference = "inverse of numDeriv Hessian of l(rho, log c, log phi)"
  )
}

information <- timed("information checks", {
  rbind(
    check_population_information(2L, c(0.7, 0.3)),
    check_population_information(3L, c(0.6, 0.3, 0.1)),
    do.call(rbind, lapply(2:4, check_single_cell_information))
  )
})

# --- 2. Coverage of the interior intervals -------------------------------

population_coverage <- function(n_celltypes, p, n_columns = 5L) {
  moments <- population_moments(n_celltypes)
  mu <- moments$mean_profiles
  sigma <- moments$covariance_matrices
  names(p) <- colnames(mu)
  methods <- c("ilr", "wald", "firth", "profile")
  rows <- lapply(seq_len(n_rep), function(r) {
    y <- withr::with_seed(seed + r, {
      simulate_bulk_mixture(mu, sigma, p = p, n = n_columns)$Y
    })
    do.call(
      rbind,
      lapply(methods, function(method) {
        fit <- tryCatch(
          suppressWarnings(confint_decovart(
            y,
            mu,
            sigma,
            method = method,
            level = level,
            seed = seed + r
          )),
          error = function(e) NULL
        )
        interval <- if (is.null(fit)) {
          matrix(NA_real_, length(p), 2L)
        } else {
          fit$interval
        }
        data.frame(
          model = "population",
          scenario = sprintf("J = %d, N = %d", n_celltypes, n_columns),
          method = method,
          celltype = names(p),
          truth = p,
          covered = covers(interval, p),
          width = interval[, 2L] - interval[, 1L]
        )
      })
    )
  })
  do.call(rbind, rows)
}

single_cell_coverage <- function(scenario_name) {
  design <- single_cell_scenario_design()
  scenario <- design$scenarios[[scenario_name]]
  ref <- single_cell_reference(2L, scenario$cell_size)
  mu <- ref$mean_signature_matrix
  sigma <- ref$Sigma
  sims <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = scenario$proportions,
    n_cells = design$n_cells,
    platform_scale = 1,
    n_samples = n_rep,
    seed = seed
  )
  methods <- expand.grid(
    method = c("ilr", "wald", "profile"),
    variance = c("ml", "reml"),
    stringsAsFactors = FALSE
  )
  rows <- lapply(seq_len(n_rep), function(r) {
    y <- sims$Y[, r]
    truth <- sims$p_realised[, r]
    scale_truth <- c(sims$mean_scale[[r]], sims$variance_scale[[r]])
    do.call(
      rbind,
      lapply(seq_len(nrow(methods)), function(k) {
        method <- methods$method[[k]]
        variance <- methods$variance[[k]]
        fit <- tryCatch(
          suppressWarnings(confint_single_cell(
            y,
            mu,
            sigma,
            method = method,
            level = level,
            variance = variance,
            seed = seed + r
          )),
          error = function(e) NULL
        )
        interval <- if (is.null(fit)) {
          matrix(NA_real_, length(truth), 2L)
        } else {
          fit$interval
        }
        scales <- if (is.null(fit) || is.null(fit$nuisance$interval)) {
          matrix(NA_real_, 2L, 2L)
        } else {
          fit$nuisance$interval
        }
        data.frame(
          model = "single-cell",
          scenario = scenario_name,
          method = paste(method, variance),
          celltype = c(names(truth), "mean_scale c", "variance_scale phi"),
          truth = c(truth, scale_truth),
          covered = c(covers(interval, truth), covers(scales, scale_truth)),
          width = c(
            interval[, 2L] - interval[, 1L],
            scales[, 2L] - scales[, 1L]
          )
        )
      })
    )
  })
  do.call(rbind, rows)
}

coverage <- rbind(
  timed("population J = 2", population_coverage(2L, c(0.7, 0.3))),
  timed("population J = 3", population_coverage(3L, c(0.6, 0.3, 0.1))),
  do.call(
    rbind,
    lapply(names(single_cell_scenario_design()$scenarios), function(nm) {
      timed(paste("single-cell", nm), single_cell_coverage(nm))
    })
  )
)

# --- 3. Bootstrap spot checks --------------------------------------------

bootstrap_spot_check <- function() {
  moments <- population_moments(2L)
  mu <- moments$mean_profiles
  sigma <- moments$covariance_matrices
  p <- stats::setNames(c(0.7, 0.3), colnames(mu))
  n_cells_ref <- 60L
  rows <- lapply(seq_len(n_rep_boot), function(r) {
    withr::with_seed(seed + r, {
      y <- simulate_bulk_mixture(mu, sigma, p = p, n = 5L)$Y
      parametric <- confint_decovart(
        y,
        mu,
        sigma,
        method = "bootstrap",
        n_boot = 100L
      )
      # Reference profiles drawn from the true moments: the reference
      # bootstrap also propagates their sampling error.
      refs <- lapply(seq_along(p), function(j) {
        t(MASS::mvrnorm(n_cells_ref, mu[, j], sigma[,, j]))
      })
      names(refs) <- names(p)
      reference <- confint_decovart(
        y,
        method = "reference",
        reference_profiles = refs,
        n_boot = 100L
      )
      data.frame(
        model = "population",
        scenario = "J = 2, N = 5",
        method = c("bootstrap", "bootstrap", "reference", "reference"),
        celltype = rep(names(p), 2L),
        truth = rep(p, 2L),
        covered = c(
          covers(parametric$interval, p),
          covers(reference$interval, p)
        ),
        width = c(
          parametric$interval[, 2L] - parametric$interval[, 1L],
          reference$interval[, 2L] - reference$interval[, 1L]
        )
      )
    })
  })
  do.call(rbind, rows)
}

bootstrap <- timed("bootstrap spot checks", bootstrap_spot_check())

# --- 4. Size of boundary likelihood-ratio tests --------------------------

boundary_size <- function() {
  cases <- list(
    list(J = 2L, p = c(1, 0), null = c(celltype_2 = 0)),
    list(J = 3L, p = c(0.6, 0.4, 0), null = c(celltype_3 = 0)),
    list(J = 3L, p = c(1, 0, 0), null = c(celltype_1 = 1))
  )
  rows <- lapply(cases, function(case) {
    moments <- population_moments(case$J)
    mu <- moments$mean_profiles
    sigma <- moments$covariance_matrices
    tests <- lapply(seq_len(n_rep), function(r) {
      y <- withr::with_seed(seed + r, {
        simulate_bulk_mixture(mu, sigma, p = case$p, n = 5L)$Y
      })
      suppressWarnings(lrt_decovart(y, mu, sigma, null_value = case$null))
    })
    p_values <- vapply(tests, function(t) t$p_value, numeric(1))
    data.frame(
      model = "population",
      null = sprintf(
        "J = %d, %s = %g",
        case$J,
        names(case$null),
        case$null
      ),
      n_boundary = tests[[1L]]$n_boundary,
      rejection_rate = mean(p_values < 1 - level),
      mc_se = mc_se(mean(p_values < 1 - level), length(p_values))
    )
  })
  design <- single_cell_scenario_design()
  ref <- single_cell_reference(2L, c(1, 5))
  sims <- simulate_bulk_from_cells(
    ref$generative_mean,
    ref$generative_sigma,
    p = c(1, 0),
    n_cells = design$n_cells,
    platform_scale = 1,
    n_samples = n_rep,
    seed = seed
  )
  sc_p <- vapply(
    seq_len(n_rep),
    function(r) {
      suppressWarnings(lrt_single_cell(
        sims$Y[, r],
        ref$mean_signature_matrix,
        ref$Sigma,
        null_value = c(celltype_2 = 0)
      ))$p_value
    },
    numeric(1)
  )
  rbind(
    do.call(rbind, rows),
    data.frame(
      model = "single-cell",
      null = "J = 2, celltype_2 = 0",
      n_boundary = 1L,
      rejection_rate = mean(sc_p < 1 - level),
      mc_se = mc_se(mean(sc_p < 1 - level), n_rep)
    )
  )
}

boundary <- timed("boundary test size", boundary_size())

# --- Summary --------------------------------------------------------------

summarise_coverage <- function(x) {
  groups <- split(
    x,
    list(x$model, x$scenario, x$method, x$celltype),
    drop = TRUE
  )
  out <- do.call(
    rbind,
    lapply(groups, function(g) {
      ok <- g$covered[!is.na(g$covered)]
      data.frame(
        model = g$model[[1L]],
        scenario = g$scenario[[1L]],
        method = g$method[[1L]],
        target = g$celltype[[1L]],
        coverage = mean(ok),
        mc_se = mc_se(mean(ok), length(ok)),
        median_width = stats::median(g$width, na.rm = TRUE),
        n_failed = sum(is.na(g$covered))
      )
    })
  )
  rownames(out) <- NULL
  out[order(out$model, out$scenario, out$method, out$target), ]
}

results <- list(
  information = information,
  coverage = summarise_coverage(coverage),
  bootstrap = summarise_coverage(bootstrap),
  boundary = boundary,
  n_rep = n_rep,
  n_rep_boot = n_rep_boot,
  level = level,
  seed = seed,
  session = sessioninfo::session_info()
)
saveRDS(results, file.path(output_dir, "ci_checks.rds"))

options(width = 120)
print(results$information, digits = 3, row.names = FALSE)
print(results$coverage, digits = 3, row.names = FALSE)
print(results$bootstrap, digits = 3, row.names = FALSE)
print(results$boundary, digits = 3, row.names = FALSE)
