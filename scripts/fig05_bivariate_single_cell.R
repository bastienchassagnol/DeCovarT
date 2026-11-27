###############################################################################
###                                                                         ###
###     FIGURE 05 - BIVARIATE TOY, SINGLE-CELL MODEL (article Section 2.4)  ###
###       J = 2 cell types, G = 2 genes, fixed counts, known scales         ###
###                                                                         ###
###############################################################################
#
# Run from the repository root; stdout and stderr go to logs/ (see
# scripts/run_figure_benchmarks.sh for the OpenBLAS settings). Workers load
# the installed DeCovarT, so run devtools::install() first:
#
#   bash scripts/run_figure_benchmarks.sh fig05
#
# The single-cell twin of scripts/fig02_bivariate_toy.R, built so that the
# covariance law is the only difference between the two. A bulk sums exactly
# m_j = n p_j cells of type j (no multinomial draw), each N(mu_j, Sigma_W,j),
# with Sigma_W,j = n R_j. The bulk per cell z = y / n is then
#   z ~ N(mu p, p_1 R_1 + p_2 R_2)            (single cell, linear in p)
# against fig02's
#   z ~ N(mu p, p_1^2 R_1 + p_2^2 R_2)        (population, quadratic in p),
# with the same means, compositions and unit-diagonal R_j. n cancels from
# the likelihood of z; the scales c = phi = n are known and given to every
# scale-dependent method, and only p_1 is estimated. The target is the
# design composition p.
#
# Factorial design (486 scenarios)
#  Centroid gap            small: mu = (20, 22) / (22, 20);
#                          large: mu = (20, 40) / (40, 20)  (as fig02)
#  Proportions             balanced (1/2, 1/2); moderately unbalanced
#                          (17/20, 3/20); highly unbalanced (99/100, 1/100)
#  Correlation, type 1     -0.8 to 0.8 by 0.2 (9 levels)
#  Correlation, type 2     -0.8 to 0.8 by 0.2 (9 levels)
#  Fixed                   n = 1,000 cells (integer counts for every p),
#                          Sigma_W,j = n R_j, platform scale a = 1
#  Scenario ID             SB{index}_{Ba|Mo|Hi}_{Sm|Lg}
#  Boundary null           8 extra scenarios at p = (1, 0):
#                          {Sm, Lg} x 4 corner correlation pairs
#
# Methods (same bulks for every method)
#  lsei                  DeconRNASeq on z = y / c (simplex-constrained LS)
#  NNLS                  non-negative least squares, free scale
#  Newton-Raphson        DeCovarT, trust-region Newton in ILR coordinates
#  L-BFGS-B              DeCovarT, L-BFGS-B on the closure of the orthant
#  Marquardt-Levenberg   DeCovarT, marqLevAlg in ILR coordinates
#  Firth                 DeCovarT with the Jeffreys penalty
#  Scalar MLE            reference fit: stats::optimize() on p_1 in (0, 1)
#                        plus both endpoints, the highest likelihood kept
#
# Inference
#  Wald, ILR, oracle      every DeCovarT solver (Wald and ILR as
#                         diagnostics); likelihood-ratio inversion at the
#                         known scales (`profile`) under Newton-Raphson
#  Parametric bootstrap   Newton-Raphson, B = 999, on 12 scenarios: small
#                         gap x 3 compositions x the corner pairs (0, 0),
#                         (-0.8, -0.8), (0.8, 0.8), (-0.8, 0.8)
#  Boundary LRT           p_2 = 0 at the known scales: chi-bar-square
#                         calibration on every bulk (power); size on the
#                         p = (1, 0) scenarios, chi-bar-square against a
#                         Monte Carlo test under the (simple) null
#  Paired population fit  the same standard-normal draws mapped through
#                         fig02's covariance p_1^2 R_1 + p_2^2 R_2, fitted
#                         by the scalar MLE of that law and by lsei
#
# Outputs (output/fig05/)
#  fig05_benchmark.rds   settings, design, truth, every metric table, the
#                        paired population fits, the boundary LRTs
#  fig05_plot_data.rds   the data behind every figure page, keyed by stem
#  density_visualisations/, performance_visualisations/: one SVG and one
#  RDS (its plot data) per page
#
# Sourcing this file from a test or a vignette defines the builders only.
###############################################################################

FIG05_N_CELLS <- 1000L

corr2 <- function(rho) matrix(c(1, rho, rho, 1), 2L, 2L)

#' Build the bivariate single-cell scenario configuration
#'
#' @return Tibble with one row per scenario, its descriptors and the
#'   `true_theta` list column read by [run_simulation_benchmark()].
build_bivariate_single_cell_config <- function(
  proportions = list(
    "balanced" = c(0.5, 0.5),
    "moderately unbalanced" = c(0.85, 0.15),
    "highly unbalanced" = c(0.99, 0.01)
  ),
  signature_matrices = list(
    "small_CLD" = matrix(c(20, 22, 22, 20), nrow = 2),
    "large_CLD" = matrix(c(20, 40, 40, 20), nrow = 2)
  ),
  corr_sequence = seq(-0.8, 0.8, 0.2),
  n_cells = FIG05_N_CELLS
) {
  proportion_list <- proportions
  genes <- c("gene_1", "gene_2")
  types <- c("celltype_1", "celltype_2")
  signature_matrices <- lapply(signature_matrices, function(mu) {
    dimnames(mu) <- list(genes, types)
    mu
  })
  design <- tidyr::expand_grid(
    centroids = names(signature_matrices),
    proportions = names(proportion_list),
    correlation_celltype1 = corr_sequence,
    correlation_celltype2 = corr_sequence
  )
  p_code <- c(
    "balanced" = "Ba",
    "moderately unbalanced" = "Mo",
    "highly unbalanced" = "Hi",
    "absent" = "Ab"
  )
  purrr::pmap(
    design,
    function(
      centroids,
      proportions,
      correlation_celltype1,
      correlation_celltype2
    ) {
      p <- stats::setNames(proportion_list[[proportions]], types)
      mu <- signature_matrices[[centroids]]
      R <- array(0, c(2L, 2L, 2L), dimnames = list(genes, genes, types))
      R[,, 1L] <- corr2(correlation_celltype1)
      R[,, 2L] <- corr2(correlation_celltype2)
      d <- mu[, 1L] - mu[, 2L]
      fisher_sc <- fisher_split(p, mu, R, "linear")
      fisher_pop <- fisher_split(p, mu, R, "quadratic")
      # Overlap of the component laws N(mu_j, R_j): the same as fig02's.
      overlap <- tryCatch(
        signif(MixSim::overlap(Pi = p, Mu = t(mu), S = R)$BarOmega, 3),
        error = function(e) NA_real_
      )
      tibble::tibble(
        centroids = centroids,
        correlation_celltype1 = correlation_celltype1,
        correlation_celltype2 = correlation_celltype2,
        overlap = overlap,
        covariance_contrast = abs(correlation_celltype1 - correlation_celltype2),
        var_parallel_1 = sum(d * (R[,, 1L] %*% d)) / sum(d^2),
        var_parallel_2 = sum(d * (R[,, 2L] %*% d)) / sum(d^2),
        entropy = round(compute_shannon_entropy(p), 3),
        proportions = proportions,
        fisher_mean_sc = fisher_sc[["mean"]],
        fisher_cov_sc = fisher_sc[["cov"]],
        fisher_mean_pop = fisher_pop[["mean"]],
        fisher_cov_pop = fisher_pop[["cov"]],
        true_theta = list(list(
          p = p,
          mu = mu,
          sigma = n_cells * R,
          n_cells = n_cells,
          platform_scale = 1
        ))
      )
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      ID = sprintf(
        "SB%d_%s_%s",
        dplyr::row_number(),
        p_code[.data$proportions],
        c(small_CLD = "Sm", large_CLD = "Lg")[.data$centroids]
      ),
      .before = 1L
    )
}

# Mixture covariance on the per-cell scale and its derivatives in p_1:
# linear (single cell) V = p_1 R_1 + p_2 R_2, quadratic (population)
# V = p_1^2 R_1 + p_2^2 R_2.
mixture_covariance <- function(p1, R, law = c("linear", "quadratic")) {
  law <- match.arg(law)
  R1 <- R[,, 1L]
  R2 <- R[,, 2L]
  switch(
    law,
    linear = list(
      V = p1 * R1 + (1 - p1) * R2,
      dV = R1 - R2,
      d2V = 0 * R1
    ),
    quadratic = list(
      V = p1^2 * R1 + (1 - p1)^2 * R2,
      dV = 2 * p1 * R1 - 2 * (1 - p1) * R2,
      d2V = 2 * (R1 + R2)
    )
  )
}

#' Fisher information of p_1, split into its mean and covariance parts
#'
#' For z ~ N(mu p, V(p)), d = mu_1 - mu_2:
#' I = d' V^-1 d + tr[(V^-1 V')^2] / 2. Under the linear law V' = R_1 - R_2
#' does not depend on p, and is zero when R_1 = R_2.
fisher_split <- function(p, mu, R, law = c("linear", "quadratic")) {
  cov <- mixture_covariance(p[[1L]], R, law)
  d <- mu[, 1L] - mu[, 2L]
  a <- solve(cov$V, cov$dV)
  c(mean = sum(d * solve(cov$V, d)), cov = 0.5 * sum(diag(a %*% a)))
}

#' Log-likelihood of p_1 for one bulk per cell, with analytic derivatives
#'
#' z ~ N(mu_2 + p_1 d, V(p_1)), without the constant -log(2 pi). Returns the
#' log-likelihood, the score, the observed second derivative and the two
#' parts of the Fisher information.
sc_terms <- function(p1, z, mu, R, law = c("linear", "quadratic")) {
  stopifnot(length(p1) == 1L, p1 >= 0, p1 <= 1)
  cov <- mixture_covariance(p1, R, law)
  d <- mu[, 1L] - mu[, 2L]
  r <- z - (mu[, 2L] + p1 * d)
  U <- chol(cov$V)
  A <- chol2inv(U)
  Ar <- drop(A %*% r)
  Ad <- drop(A %*% d)
  B <- A %*% cov$dV
  BAr <- drop(B %*% Ar)
  list(
    loglik = -sum(log(diag(U))) - 0.5 * sum(r * Ar),
    score = sum(d * Ar) + 0.5 * sum(Ar * drop(cov$dV %*% Ar)) -
      0.5 * sum(diag(B)),
    hessian = -sum(d * Ad) -
      2 * sum(Ad * drop(cov$dV %*% Ar)) -
      sum(BAr * drop(cov$dV %*% Ar)) +
      0.5 * sum(Ar * drop(cov$d2V %*% Ar)) +
      0.5 * sum(diag(B %*% B)) -
      0.5 * sum(A * cov$d2V),
    fisher_mean = sum(d * Ad),
    fisher_cov = 0.5 * sum(diag(B %*% B))
  )
}

#' Reference maximum-likelihood fit of p_1 on [0, 1]
#'
#' `stats::optimize()` on the interior, then both endpoints: the realised
#' likelihood need not be concave, and the maximum may sit on a face.
fit_sc_scalar <- function(z, mu, R, law = "linear", tol = 1e-10) {
  ll <- function(p1) sc_terms(p1, z, mu, R, law)$loglik
  interior <- stats::optimize(ll, c(0, 1), maximum = TRUE, tol = tol)
  candidates <- c(0, interior$maximum, 1)
  values <- vapply(candidates, ll, numeric(1))
  k <- which.max(values)
  list(
    p = stats::setNames(c(candidates[[k]], 1 - candidates[[k]]), colnames(mu)),
    loglik = values[[k]]
  )
}

# The scalar MLE as a benchmark method: y ~ N(c mu p, phi S_W(p)) is
# z = y / c ~ N(mu p, (phi / c^2) S_W(p)).
scalar_mle_plug_in <- function(
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  fit_sc_scalar(
    y / mean_scale,
    mean_signature_matrix,
    Sigma * variance_scale / mean_scale^2
  )$p
}

# The scalar MLE with its helpers in its own environment: future workers
# receive the whole closure, as they do not see this script's functions.
scalar_mle_method <- function() {
  env <- new.env(parent = globalenv())
  for (name in c(
    "mixture_covariance",
    "sc_terms",
    "fit_sc_scalar",
    "scalar_mle_plug_in"
  )) {
    fun <- get(name)
    environment(fun) <- env
    assign(name, fun, envir = env)
  }
  env$scalar_mle_plug_in
}

# Closed-form least squares on the segment between the two signatures.
fit_lsei_scalar <- function(z, mu) {
  d <- mu[, 1L] - mu[, 2L]
  p1 <- min(1, max(0, sum(d * (z - mu[, 2L])) / sum(d^2)))
  stats::setNames(c(p1, 1 - p1), colnames(mu))
}

# lsei reads the bulk per cell, y / c.
lsei_plug_in <- function(y, mean_signature_matrix, mean_scale) {
  deconvolute_ratios_deconrnaseq(y / mean_scale, mean_signature_matrix)
}

#' Single-cell likelihood with fixed cell counts and known scales
#'
#' A copy of `single_cell_likelihood(n, n)` whose generator sums exactly
#' m_j = n p_j cells: y ~ N(n mu p, n S_W(p)), the target is the design p
#' and the scales are c = phi = n in every bulk. The multinomial generator
#' of `simulate_bulk_from_cells()` (used by fig06) is not called.
fixed_count_single_cell_likelihood <- function(n_cells = FIG05_N_CELLS) {
  likelihood <- single_cell_likelihood(n_cells, n_cells)
  likelihood$simulate <- function(mean, Sigma, theta, n) {
    p <- as.numeric(theta$p)
    s_w <- apply(sweep(Sigma, 3L, p, "*"), c(1L, 2L), sum)
    eps <- matrix(stats::rnorm(nrow(mean) * n), nrow(mean), n)
    Y <- theta$n_cells * drop(mean %*% p) +
      t(chol(theta$n_cells * s_w)) %*% eps
    samples <- paste0("sample_", seq_len(n))
    dimnames(Y) <- list(rownames(mean), samples)
    list(
      Y = Y,
      truth = matrix(p, length(p), n, dimnames = list(colnames(mean), samples)),
      nuisance = rbind(
        mean_scale = rep(theta$n_cells, n),
        variance_scale = rep(theta$n_cells, n)
      )
    )
  }
  likelihood
}

#' Methods and intervals of the bivariate single-cell benchmark
#'
#' @param intervals If `FALSE` (tests), no confidence intervals.
#' @param bootstrap If `TRUE`, Newton-Raphson only, with the bootstrap.
bivariate_single_cell_deconvolution_functions <- function(
  itmax = 100L,
  epsilon = 1e-6,
  intervals = TRUE,
  bootstrap = FALSE,
  n_boot = 999L
) {
  solver <- list(epsilon = epsilon, itmax = itmax)
  menu <- function(methods, optimiser = "newton") {
    if (!isTRUE(intervals)) {
      return(NULL)
    }
    list(methods = methods, optimiser = optimiser)
  }
  decovart <- function(optimiser, ci_methods) {
    list(
      FUN = deconvolute_single_cell,
      additional_parameters = c(solver, optimiser = optimiser),
      scales = "true",
      ci = menu(ci_methods, optimiser)
    )
  }
  if (isTRUE(bootstrap)) {
    newton <- decovart("newton", NULL)
    newton$ci <- list(methods = "bootstrap", optimiser = "newton", n_boot = n_boot)
    return(list("Newton-Raphson" = newton))
  }
  list(
    "lsei" = list(FUN = lsei_plug_in, scales = "true"),
    "NNLS" = list(FUN = deconvolute_ratios_nnls),
    "Newton-Raphson" = decovart("newton", c("wald", "ilr", "profile", "oracle")),
    "L-BFGS-B" = decovart("lbfgsb", c("wald", "ilr")),
    "Marquardt-Levenberg" = decovart("marquardt", c("wald", "ilr")),
    "Firth" = list(
      FUN = fit_firth_single_cell,
      scales = "true",
      ci = menu(c("wald", "ilr"))
    ),
    "Scalar MLE" = list(FUN = scalar_mle_method(), scales = "true")
  )
}

# Per-cell bulk z and per-cell covariances R of one scenario.
.per_cell <- function(theta, bulk) {
  list(
    z = bulk / theta$n_cells,
    mu = theta$mu,
    R = theta$sigma / theta$n_cells
  )
}

# One list entry per scenario (ID, design p, per-cell z, mu, R). Workers get
# these slices only: a closure over the benchmark would ship every scenario's
# bulks to every worker.
.scenario_inputs <- function(benchmark, config) {
  lapply(seq_len(nrow(config)), function(i) {
    theta <- config$true_theta[[i]]
    c(
      list(ID = config$ID[[i]], p = as.numeric(theta$p)),
      .per_cell(theta, benchmark$data[[i]]$bulk)
    )
  })
}

#' Paired population fits on common random numbers
#'
#' Each saved bulk is z = mu p + L_sc eps with L_sc the Cholesky factor of
#' p_1 R_1 + p_2 R_2; eps is recovered exactly and mapped through fig02's
#' covariance p_1^2 R_1 + p_2^2 R_2. Both bulks are fitted by the scalar MLE
#' of their own law and by lsei, so the difference between the two models
#' carries little Monte Carlo noise.
paired_one <- function(input) {
  p <- input$p
  m <- drop(input$mu %*% p)
  l_sc <- t(chol(mixture_covariance(p[[1L]], input$R, "linear")$V))
  l_pop <- t(chol(mixture_covariance(p[[1L]], input$R, "quadratic")$V))
  eps <- forwardsolve(l_sc, input$z - m)
  z_pop <- m + l_pop %*% eps
  fits <- lapply(seq_len(ncol(input$z)), function(s) {
    c(
      sc_mle = fit_sc_scalar(input$z[, s], input$mu, input$R, "linear")$p[[1L]],
      sc_lsei = fit_lsei_scalar(input$z[, s], input$mu)[[1L]],
      pop_mle = fit_sc_scalar(z_pop[, s], input$mu, input$R, "quadratic")$p[[
        1L
      ]],
      pop_lsei = fit_lsei_scalar(z_pop[, s], input$mu)[[1L]]
    )
  })
  fits <- do.call(rbind, fits)
  tibble::tibble(
    ID = input$ID,
    sample_id = rep(colnames(input$z), 2L),
    p1_true = p[[1L]],
    method = rep(c("Scalar MLE", "lsei"), each = nrow(fits)),
    p1_single_cell = c(fits[, "sc_mle"], fits[, "sc_lsei"]),
    p1_population = c(fits[, "pop_mle"], fits[, "pop_lsei"])
  )
}

paired_population_table <- function(benchmark, config, cores = 1L) {
  .map_scenarios(.scenario_inputs(benchmark, config), paired_one, cores)
}

.map_scenarios <- function(inputs, fun, cores) {
  out <- if (cores > 1L) {
    furrr::future_map(inputs, fun, .options = furrr::furrr_options(seed = TRUE))
  } else {
    lapply(inputs, fun)
  }
  dplyr::bind_rows(out)
}

# LR statistic of p_2 = 0 (p_1 = 1) at known scales.
lrt_scalar <- function(z, mu, R) {
  full <- fit_sc_scalar(z, mu, R)$loglik
  max(0, 2 * (full - sc_terms(1, z, mu, R)$loglik))
}

#' Null scenarios of the boundary test: p = (1, 0)
#'
#' Both centroid gaps crossed with the four corner correlation pairs.
boundary_null_config <- function(
  corners = list(c(0, 0), c(-0.8, -0.8), c(0.8, 0.8), c(-0.8, 0.8))
) {
  config <- build_bivariate_single_cell_config(
    proportions = list("absent" = c(1, 0)),
    corr_sequence = c(-0.8, 0, 0.8)
  )
  keep <- paste(config$correlation_celltype1, config$correlation_celltype2) %in%
    vapply(corners, paste, character(1), collapse = " ")
  config <- config[keep, ]
  config$ID <- sprintf(
    "SB0%d_Ab_%s",
    seq_len(nrow(config)),
    c(small_CLD = "Sm", large_CLD = "Lg")[config$centroids]
  )
  config
}

#' Boundary likelihood-ratio tests of p_2 = 0
#'
#' `power`: every bulk of the main run, chi-bar-square p-value
#' (1/2 chi2_0 + 1/2 chi2_1). `size`: `n_null` bulks of each p = (1, 0)
#' scenario, with the chi-bar-square p-value and the Monte Carlo p-value
#' (1 + #{T_b >= T}) / (B + 1) against `b_null` statistics simulated under
#' the null. The null is simple (no nuisance parameter), so one null sample
#' per scenario serves every bulk; the Monte Carlo test is then exact up to
#' its Monte Carlo error. `check`: package `lrt_single_cell()` at the known
#' scales on the first `n_check` null bulks, against the scalar statistic.
.chibar <- function(t) chi_bar_square_pvalue(t, n_boundary = 1L)

boundary_power_one <- function(input) {
  statistic <- apply(input$z, 2L, lrt_scalar, mu = input$mu, R = input$R)
  tibble::tibble(
    ID = input$ID,
    sample_id = colnames(input$z),
    statistic = statistic,
    p_value = vapply(statistic, .chibar, numeric(1))
  )
}

boundary_size_one <- function(input) {
  theta <- input$theta
  draw <- function(n) {
    fixed_count_single_cell_likelihood(theta$n_cells)$simulate(
      theta$mu,
      theta$sigma,
      theta,
      n
    )$Y
  }
  withr::with_seed(input$seed, {
    observed <- draw(input$n_null)
    reference <- draw(input$b_null)
  })
  R <- theta$sigma / theta$n_cells
  stat <- function(Y) {
    apply(Y / theta$n_cells, 2L, lrt_scalar, mu = theta$mu, R = R)
  }
  t_obs <- stat(observed)
  t_ref <- stat(reference)
  check <- vapply(
    seq_len(min(input$n_check, input$n_null)),
    function(s) {
      suppressWarnings(lrt_single_cell(
        observed[, s],
        theta$mu,
        theta$sigma,
        null_value = c(celltype_2 = 0),
        mean_scale = theta$n_cells,
        variance_scale = theta$n_cells
      ))$statistic
    },
    numeric(1)
  )
  tibble::tibble(
    ID = input$ID,
    sample_id = seq_len(input$n_null),
    statistic = t_obs,
    p_chibar = vapply(t_obs, .chibar, numeric(1)),
    p_monte_carlo = (1 + vapply(t_obs, function(t) sum(t_ref >= t), 0)) /
      (input$b_null + 1),
    package_statistic = c(check, rep(NA_real_, input$n_null - length(check)))
  )
}

boundary_lrt_tables <- function(
  benchmark,
  config,
  null_config,
  n_null = 5000L,
  b_null = 4999L,
  n_check = 200L,
  seed = 1L,
  cores = 1L
) {
  power <- .map_scenarios(
    .scenario_inputs(benchmark, config),
    boundary_power_one,
    cores
  )
  size_inputs <- lapply(seq_len(nrow(null_config)), function(i) {
    list(
      ID = null_config$ID[[i]],
      theta = null_config$true_theta[[i]],
      seed = seed + i,
      n_null = n_null,
      b_null = b_null,
      n_check = n_check
    )
  })
  size <- .map_scenarios(size_inputs, boundary_size_one, cores)
  list(
    power = power,
    size = size,
    null_config = null_config[setdiff(names(null_config), "true_theta")],
    settings = list(n_null = n_null, b_null = b_null, seed = seed)
  )
}

if (
  interactive() ||
    any(grepl(
      "fig05_bivariate_single_cell\\.R$",
      sub(
        "^--file=",
        "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      )
    ))
) {
  # Run size and output. FIG05_POSTPROCESS_ONLY=1 redraws the figures from
  # the saved benchmark.
  N_REPLICATES <- as.integer(Sys.getenv("FIG05_N_REPLICATES", "2000"))
  N_BOOT <- 999L
  N_NULL <- 5000L
  B_NULL <- 4999L
  # Half of the machine's cores (44 of 88).
  CORES <- parallel::detectCores() %/% 2L
  POSTPROCESS_ONLY <- identical(Sys.getenv("FIG05_POSTPROCESS_ONLY"), "1")
  OUT_DIR <- file.path("output", "fig05")
  BENCHMARK_FILE <- file.path(OUT_DIR, "fig05_benchmark.rds")
  CHECKPOINT_FILE <- file.path(OUT_DIR, "fig05_checkpoint_main.rds")
  RESUME <- identical(Sys.getenv("FIG05_RESUME"), "1")
  SEED <- 20261125L
  ITMAX <- 100L
  EPSILON <- 1e-6

  library(DeCovarT)
  local({
    a <- 0.5^abs(outer(seq_len(64L), seq_len(64L), "-"))
    if (
      max(abs(crossprod(chol(a)) - a)) > 1e-10 ||
        max(abs(solve(a) %*% a - diag(64L))) > 1e-8
    ) {
      stop(
        "BLAS returns wrong chol()/solve(); set OPENBLAS_CORETYPE=Haswell ",
        "before starting R.",
        call. = FALSE
      )
    }
  })
  source(file.path("scripts", "benchmark_metrics", "benchmark_artefacts.R"))
  source(file.path("scripts", "visualisation", "load_visualisation.R"))
  stopifnot(requireNamespace("MixSim", quietly = TRUE))

  cli::cli_h1("Figure 05: bivariate toy, single-cell model")

  if (!POSTPROCESS_ONLY) {
    config <- build_bivariate_single_cell_config()
    cli::cli_alert_info(
      "{.val {nrow(config)}} scenarios, N = {.val {N_REPLICATES}}, {.val {CORES}} workers."
    )
    if (CORES > 1L) {
      future::plan(future::multisession, workers = CORES)
    }
    started <- proc.time()[["elapsed"]]
    # Main run and bootstrap are checkpointed: FIG05_RESUME=1 restarts from
    # the paired and boundary stages.
    if (RESUME && file.exists(CHECKPOINT_FILE)) {
      out <- readRDS(CHECKPOINT_FILE)
    } else {
    out <- run_simulation_benchmark(
      scenario_config = config,
      deconvolution_functions = bivariate_single_cell_deconvolution_functions(
        itmax = ITMAX,
        epsilon = EPSILON
      ),
      n = N_REPLICATES,
      cores = CORES,
      verbose = TRUE,
      likelihood = fixed_count_single_cell_likelihood(),
      return_data = TRUE,
      # Scenario i runs under seed SEED + i: the bootstrap pass below
      # redraws the same bulks.
      seed = SEED
    )

    boot_ids <- which(
      config$centroids == "small_CLD" &
        paste(config$correlation_celltype1, config$correlation_celltype2) %in%
          c("0 0", "-0.8 -0.8", "0.8 0.8", "-0.8 0.8")
    )
    stopifnot(length(boot_ids) == 12L)
    cli::cli_alert_info("Parametric bootstrap on {length(boot_ids)} scenarios.")
    boot <- run_simulation_benchmark(
      scenario_config = config,
      deconvolution_functions = bivariate_single_cell_deconvolution_functions(
        itmax = ITMAX,
        epsilon = EPSILON,
        bootstrap = TRUE,
        n_boot = N_BOOT
      ),
      n = N_REPLICATES,
      cores = CORES,
      verbose = TRUE,
      likelihood = fixed_count_single_cell_likelihood(),
      seed = SEED,
      scenarios = boot_ids
    )
    out$intervals <- dplyr::bind_rows(out$intervals, boot$intervals)
    rm(boot)
    saveRDS(out, CHECKPOINT_FILE)
    }

    cli::cli_alert_info("Paired population fits on the same draws.")
    out$paired <- paired_population_table(out, config, cores = CORES)

    cli::cli_alert_info("Boundary likelihood-ratio tests of p_2 = 0.")
    out$boundary_lrt <- boundary_lrt_tables(
      out,
      config,
      boundary_null_config(),
      n_null = N_NULL,
      b_null = B_NULL,
      seed = SEED + 10000L,
      cores = CORES
    )

    cli::cli_alert_success(
      "Benchmark finished in {.val {round(proc.time()[['elapsed']] - started)}} s."
    )
    save_benchmark(
      out,
      config,
      benchmark_meta(SEED, N_REPLICATES),
      BENCHMARK_FILE
    )
  }

  benchmark <- readRDS(BENCHMARK_FILE)
  plot_data <- figure_plot_data(benchmark, figure = "fig05")
  saveRDS(plot_data, file.path(OUT_DIR, "fig05_plot_data.rds"))
  render_figure_books(plot_data, OUT_DIR)
  # Pages shown in vignettes/R3-bivariate-single-cell.qmd (pkgdown only).
  shown <- c(
    "density_visualisations/purified_density/small_centroid_gap_balanced",
    "density_visualisations/loglik_mean_shift_slice/small_centroid_gap_moderately_unbalanced",
    "density_visualisations/loglik_correlation_slice/small_centroid_gap_balanced",
    "performance_visualisations/heatmap_rmse_ratio_lsei/small_centroid_gap_highly_unbalanced",
    "performance_visualisations/rmse_overlap_scatter/small_centroid_gap",
    "performance_visualisations/heatmap_rmse_ratio_population/small_centroid_gap_moderately_unbalanced",
    "performance_visualisations/boundary_lrt/boundary_lrt"
  )
  target <- file.path("vignettes", "figures", "fig05")
  unlink(list.files(target, full.names = TRUE))
  dir.create(target, showWarnings = FALSE, recursive = TRUE)
  file.copy(
    file.path(OUT_DIR, paste0(shown, ".svg")),
    file.path(target, paste0(gsub("/", "__", shown, fixed = TRUE), ".svg")),
    overwrite = TRUE
  )
  cli::cli_alert_success("Figure 05 written to {.file {OUT_DIR}}.")
}
