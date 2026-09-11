###############################################################################
###############################################################################
###                                                                         ###
###           FIGURE 02 – BIVARIATE TOY MODEL (article 2.1)               ###
###       J = 2 cell types · G = 2 genes · full factorial design            ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Background launch from the repository root (vanilla Rscript; no CLI
# parser — hyperparameters are hard-coded below). stdout and stderr go
# to logs/:
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore scripts/fig02_bivariate_toy.R \
#     > "logs/fig02_$(date +%F)_bivariate_toy.log" 2>&1 &
#
# Redraw figures from a finished ADEMP RDS (does not refit the 972
# scenarios). Expected-Fisher Wald SEs at the true composition are
# attached to `monte_carlo` without refitting. PDFs:
#   output/fig02/density_visualisations/ and
#   output/fig02/performance_visualisations/.
#
#   FIG02_POSTPROCESS_ONLY=1 Rscript --no-save --no-restore \
#     scripts/fig02_bivariate_toy.R
#
# Article:  DeCovarT – Section 2.1 "Toy Model with two genes and two cell
#           populations" (Numerical Simulation Study).
# Vignette: vignettes/fig02-bivariate-toy.qmd
#
# Scenario builders live at the top of this file. Sourcing it from tests
# or a vignette defines those functions without running the grid:
# non-interactive `source()` does not see `--file=...fig02_bivariate_toy.R`.
#
# ── Factorial design ────────────────────────────────────────────────────────
#  Factor                  Levels
#  ─────────────────────── ───────────────────────────────────────────────────
#  Proportions             balanced (1/2, 1/2);
#                          moderately unbalanced (17/20, 3/20);
#                          highly unbalanced (99/100, 1/100)
#  Mean distance (CLD)     small: μ=(20,22)/(22,20); large: μ=(20,40)/(40,20)
#  Gene–gene corr CT1 (ρ)  −0.8 to +0.8, step 0.2  (9 levels)
#  Gene–gene corr CT2 (ρ)  −0.8 to +0.8, step 0.2  (9 levels)
#  Variance structure      homoscedastic (σ²=1,1); heteroscedastic (σ²=1,2)
#  Algorithms              NNLS, DeconRNASeq (LSEI), L-BFGS-B, gradient,
#                          Newton–Raphson, Marquardt–Levenberg, SA
#  ─────────────────────── ───────────────────────────────────────────────────
#  Total scenarios:  3 × 2 × 9 × 9 × 2 = 972
#  Replicates (n):   500
#  Scenario ID:      B{idx}_{Ho|He}_{Ba|Mo|Hi}_{Sm|Lg}
#
# ── Solver hyperparameters (pipeline; bivariate_toy_deconvolution_functions)
#  itmax      200     max. iterations (L-BFGS-B, gradient, Newton,
#                     Marquardt–Levenberg, SA)
#  epsilon    1e-4    convergence tolerance for those solvers
#  cores      1       scenario loop is sequential; do not nest workers
#  Tests use the factory defaults itmax = 10, epsilon = 1e-3.
#
# ── Usage ───────────────────────────────────────────────────────────────────
#  Rscript scripts/fig02_bivariate_toy.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/fig02/bivariate_config.rds       – slim design grid (keyed by ID)
#  output/fig02/bivariate_descriptors.rds  – scenario geometry / MixSim /
#                                            Hellinger / Jeffreys
#  output/fig02/bivariate_theta.rds        – true_theta list column
#  output/fig02/bivariate_benchmark.rds    – metrics only (regression,
#                                            monte_carlo, optimisation, call)
#  output/fig02/density_visualisations/{purified,bulk,loglik_*}.pdf
#  output/fig02/performance_visualisations/{heatmap_*,raincloud,forest,similarity,solver_dots}.pdf
#  output/fig02/ggplot_rds/*.rds            – ggplot `data` for each book
#                                            (not rgl snapshots)
###############################################################################

#' Encode a bivariate-toy scenario ID
#'
#' Pattern `B{index}_{Ho|He}_{Ba|Mo|Hi}_{Sm|Lg}`. Kept in this script so
#' tests can `source()` the builders without calling unexported helpers.
.encode_bivariate_id <- function(
  scenario_idx,
  variance,
  proportions,
  centroids
) {
  var_code <- ifelse(variance == "homoscedastic", "Ho", "He")
  p_code <- dplyr::case_when(
    proportions == "balanced" ~ "Ba",
    proportions == "moderately unbalanced" ~ "Mo",
    proportions == "highly unbalanced" ~ "Hi",
    TRUE ~ "Xx"
  )
  cld_code <- ifelse(
    grepl("small", centroids, ignore.case = TRUE),
    "Sm",
    "Lg"
  )
  paste0(
    "B",
    as.integer(scenario_idx),
    "_",
    var_code,
    "_",
    p_code,
    "_",
    cld_code
  )
}

#' Build bivariate generative-model scenario configuration
#'
#' @param proportions Named list of simplex vectors \eqn{\boldsymbol{p}}.
#' @param signature_matrices Named list of mean matrices
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{2\times 2}^{+}}.
#' @param corr_sequence Numeric sequence of within-cell-type correlations.
#' @param diagonal_terms Named list of diagonal variance templates.
#'
#' @return Tibble with `true_theta` list column plus scenario metadata
#'   for [run_simulation_benchmark()].
build_bivariate_scenario_config <- function(
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
  diagonal_terms = list(
    "homoscedastic" = c(1, 1),
    "heteroscedastic" = c(1, 2)
  )
) {
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    stop(
      "build_bivariate_scenario_config() requires MixSim.",
      call. = FALSE
    )
  }

  # Label genes and cell types on each mean signature (G x J).
  num_celltypes <- ncol(signature_matrices[[1L]])
  num_genes <- nrow(signature_matrices[[1L]])
  signature_matrices <- purrr::map(
    signature_matrices,
    function(.mean_signature_matrix) {
      dimnames(.mean_signature_matrix) <- list(
        paste0("gene_", seq_len(num_genes)),
        paste0("celltype_", seq_len(num_celltypes))
      )
      .mean_signature_matrix
    }
  )

  # Full factorial grid: centroids x p x rho_1 x rho_2 x variance.
  proportion_list <- proportions
  design <- tidyr::expand_grid(
    centroids = names(signature_matrices),
    proportion_name = names(proportion_list),
    correlation_celltype1 = corr_sequence,
    correlation_celltype2 = corr_sequence,
    variance = names(diagonal_terms)
  ) |>
    dplyr::mutate(
      scenario_idx = dplyr::row_number(),
      ID = .encode_bivariate_id(
        .data$scenario_idx,
        .data$variance,
        .data$proportion_name,
        .data$centroids
      )
    )

  purrr::pmap(
    design,
    function(
      centroids,
      proportion_name,
      correlation_celltype1,
      correlation_celltype2,
      variance,
      scenario_idx,
      ID
    ) {
      mu <- signature_matrices[[centroids]]
      p <- proportion_list[[proportion_name]]
      diag_terms <- diagonal_terms[[variance]]

      # Exchangeable correlation per cell type, then Sigma_j = D^{1/2} R D^{1/2}.
      corr_matrix <- array(
        0,
        dim = c(num_genes, num_genes, num_celltypes),
        dimnames = list(
          paste0("gene_", seq_len(num_genes)),
          paste0("gene_", seq_len(num_genes)),
          paste0("celltype_", seq_len(num_celltypes))
        )
      )
      Sigma <- corr_matrix
      corr_matrix[,, 1] <- correlation_celltype1
      corr_matrix[,, 2] <- correlation_celltype2
      for (j in seq_len(num_celltypes)) {
        diag(corr_matrix[,, j]) <- 1
        Sigma[,, j] <- sqrt(diag(diag_terms)) %*%
          corr_matrix[,, j] %*%
          sqrt(diag(diag_terms))
      }

      true_theta <- list(p = p, mu = mu, sigma = Sigma)
      overlap <- tryCatch(
        {
          overlap_fit <- MixSim::overlap(
            Pi = p,
            Mu = t(mu),
            S = Sigma
          )
          signif(overlap_fit$BarOmega, digits = 3)
        },
        error = function(e) NA_real_
      )

      tibble::tibble(
        ID = ID,
        correlation_celltype1 = correlation_celltype1,
        correlation_celltype2 = correlation_celltype2,
        overlap = overlap,
        entropy = round(compute_shannon_entropy(p), digits = 3),
        proportions = proportion_name,
        variance = variance,
        centroids = centroids,
        true_theta = list(true_theta)
      )
    }
  ) |>
    dplyr::bind_rows()
}

#' Default deconvolution solvers for the bivariate toy benchmark
#'
#' Omits CIBERSORT (too few genes for nu-SVR tuning).
#'
#' @param itmax Maximum iterations for gradient-based DeCovarT solvers.
#' @param epsilon Convergence tolerance for DeCovarT solvers.
#'
#' @return Named list suitable for [run_simulation_benchmark()].
bivariate_toy_deconvolution_functions <- function(
  itmax = 10L,
  epsilon = 1e-3
) {
  list(
    "nnls" = list(FUN = deconvolute_ratios_nnls),
    "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
    "LBFGS" = list(
      FUN = deconvolute_ratios_L_BFGS_B,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "gradient" = list(
      FUN = deconvolute_ratios_gradient_descent,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "Newton-Raphson" = list(
      FUN = deconvolute_ratios_Newton_Raphson,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "Marquardt-Levenberg" = list(
      FUN = deconvolute_ratios_Marquardt_Levenberg,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "SA" = list(
      FUN = deconvolute_ratios_simulated_annealing,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    )
  )
}

# Run the ADEMP pipeline only for an interactive Source or
# `Rscript scripts/fig02_bivariate_toy.R`. Tests copy this file into a
# withr temp directory and source it there (builders only).
if (
  interactive() ||
    any(grepl(
      "fig02_bivariate_toy\\.R$",
      sub(
        "^--file=",
        "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      )
    ))
) {
  # ==========================================================================
  # SECTION 0 · Dependencies and paths ----
  # ==========================================================================

  # Prefer the working tree over a stale user-library install.
  if (
    requireNamespace("devtools", quietly = TRUE) &&
      file.exists("DESCRIPTION")
  ) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    library(DeCovarT)
  }
  DeCovarT:::.ui_attach_script()

  stopifnot(
    requireNamespace("MixSim", quietly = TRUE),
    requireNamespace("ggplot2", quietly = TRUE)
  )

  OUT_DIR <- file.path("output", "fig02")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  .ui_h1("Figure 02 · Bivariate toy model")

  N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "500"))
  POSTPROCESS_ONLY <- identical(
    Sys.getenv("FIG02_POSTPROCESS_ONLY", "0"),
    "1"
  )

  SEED <- 20260903L
  set.seed(SEED)

  # ==========================================================================
  # SECTION 1 · GENERATIVE MODEL ----
  # ==========================================================================

  combined_path <- file.path(OUT_DIR, "bivariate_benchmark_combined.rds")
  bench_path <- file.path(OUT_DIR, "bivariate_benchmark.rds")
  config_path <- file.path(OUT_DIR, "bivariate_config.rds")

  if (isTRUE(POSTPROCESS_ONLY)) {
    .ui_info("Post-process only: skipping ADEMP refit.")
    src <- if (file.exists(combined_path)) {
      combined_path
    } else {
      bench_path
    }
    if (!file.exists(src)) {
      .ui_abort(
        "No benchmark RDS at {.file {src}}. Run the full script first."
      )
    }
    bivariate_out <- readRDS(src)
    if (file.exists(config_path)) {
      cfg_file <- readRDS(config_path)
      if ("true_theta" %in% names(cfg_file)) {
        bivariate_out$config <- cfg_file
      }
    }
    if (
      !file.exists(combined_path) &&
        !is.null(bivariate_out$optimisation)
    ) {
      .ui_info(
        "Backing up combined benchmark to {.file {combined_path}}."
      )
      file.copy(src, combined_path, overwrite = FALSE)
    }
  } else {
    .ui_info("Building scenario config (972 factorial rows).")
    scenario_config <- build_bivariate_scenario_config()
    .ui_success(
      "Config built: {.val {nrow(scenario_config)}} scenarios."
    )

    # ========================================================================
    # SECTION 2 · INFERENCE ----
    # ========================================================================

    ITMAX <- 200L
    EPSILON <- 1e-4
    deconvolution_functions <- bivariate_toy_deconvolution_functions(
      itmax = ITMAX,
      epsilon = EPSILON
    )

    .ui_info(
      "Running ADEMP benchmark with {.val {N_REPL}} replicates."
    )
    bivariate_out <- run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = deconvolution_functions,
      n = N_REPL,
      cores = 1L,
      verbose = TRUE
    )
    bivariate_out$config <- scenario_config
  }

  rewrite_ids <- isTRUE(POSTPROCESS_ONLY)
  write_simulation_artefacts(
    benchmark = bivariate_out,
    dir = OUT_DIR,
    stem = "bivariate",
    config = bivariate_out$config,
    rewrite_bivariate_id = rewrite_ids
  )
  .ui_success("Wrote split config / descriptors / theta / metrics RDS.")

  artefacts <- read_simulation_artefacts(
    OUT_DIR,
    "bivariate",
    assemble = TRUE
  )
  artefacts$theta <- readRDS(file.path(OUT_DIR, "bivariate_theta.rds"))
  artefacts <- DeCovarT:::.attach_expected_fisher_wald(artefacts)
  write_simulation_artefacts(
    benchmark = artefacts,
    dir = OUT_DIR,
    stem = "bivariate",
    config = artefacts$config
  )
  .ui_success(
    "Refreshed Wald coverage from expected Fisher at the true composition."
  )

  DENSITY_DIR <- file.path(OUT_DIR, "density_visualisations")
  PERF_DIR <- file.path(OUT_DIR, "performance_visualisations")
  GGPLOT_RDS_DIR <- file.path(OUT_DIR, "ggplot_rds")
  dir.create(DENSITY_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(PERF_DIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(GGPLOT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

  # ==========================================================================
  # SECTION 3 · VISUALISATIONS ----
  # ==========================================================================

  .ui_info("Drawing RMSE / MAE / Aitchison tile heatmaps.")
  save_bivariate_metric_heatmaps(artefacts, PERF_DIR, data_rds = GGPLOT_RDS_DIR)
  .ui_success("Saved RMSE, MAE, and Aitchison heatmap PDFs.")

  cfg <- artefacts$config
  theta_tbl <- artefacts$theta
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    .ui_info("Drawing 12-page density books (four correlation corners).")
    save_bivariate_purified_density_book(
      cfg,
      theta_tbl,
      file.path(DENSITY_DIR, "purified_density.pdf"),
      data_rds = GGPLOT_RDS_DIR
    )
    save_bivariate_bulk_density_book(
      cfg,
      theta_tbl,
      file.path(DENSITY_DIR, "bulk_density.pdf"),
      data_rds = GGPLOT_RDS_DIR
    )
    save_bivariate_loglik_surface_p_book(
      cfg,
      theta_tbl,
      file.path(DENSITY_DIR, "loglik_surface_p.pdf"),
      data_rds = GGPLOT_RDS_DIR
    )
    save_bivariate_loglik_ilr_profile_book(
      cfg,
      theta_tbl,
      file.path(DENSITY_DIR, "loglik_ilr_profile.pdf"),
      data_rds = GGPLOT_RDS_DIR
    )
    save_bivariate_loglik_rgl_book(
      cfg,
      theta_tbl,
      file.path(DENSITY_DIR, "loglik_rgl.pdf"),
      html_file = file.path(DENSITY_DIR, "loglik_rgl.html")
    )
    .ui_success("Saved density and log-likelihood PDF books.")
  } else {
    .ui_warn("{.pkg gridExtra} not available; skipping density books.")
  }

  if (requireNamespace("ggdist", quietly = TRUE)) {
    .ui_info("Drawing 12-page raincloud book (four correlation corners).")
    save_bivariate_raincloud_book(
      artefacts,
      file.path(PERF_DIR, "raincloud.pdf"),
      data_rds = GGPLOT_RDS_DIR
    )
    .ui_success("Saved {.file raincloud.pdf}.")
  } else {
    .ui_warn("{.pkg ggdist} not available; skipping raincloud book.")
  }

  .ui_info("Drawing 12-page Wald forest book (four correlation corners).")
  save_bivariate_forest_book(
    artefacts,
    file.path(PERF_DIR, "forest.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file forest.pdf} (Wald solvers only).")

  if (length(unique(artefacts$monte_carlo$algorithm)) >= 2L) {
    if (
      requireNamespace("ggdendro", quietly = TRUE) &&
        requireNamespace("cowplot", quietly = TRUE) &&
        requireNamespace("gridExtra", quietly = TRUE)
    ) {
      .ui_info("Drawing 12-page similarity book with dendrograms.")
      save_bivariate_similarity_book(
        artefacts,
        file.path(PERF_DIR, "similarity.pdf"),
        data_rds = GGPLOT_RDS_DIR
      )
      .ui_success("Saved {.file similarity.pdf}.")
    } else {
      .ui_warn(
        "Similarity book skipped (need ggdendro, cowplot, gridExtra)."
      )
    }
  }

  .ui_info("Drawing 12-page RMSE / Aitchison solver-dot book.")
  save_bivariate_solver_dots_book(
    artefacts,
    file.path(PERF_DIR, "solver_dots.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file solver_dots.pdf}.")

  .ui_success(
    "Done. Outputs in {.path {normalizePath(OUT_DIR, mustWork = FALSE)}}."
  )
}
