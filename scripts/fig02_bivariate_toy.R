###############################################################################
###############################################################################
###                                                                         ###
###           FIGURE 02 – BIVARIATE TOY MODEL (article §2.1)               ###
###       J = 2 cell types · G = 2 genes · full factorial design            ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Article:  DeCovarT – Section 2.1 "Toy Model with two genes and two cell
#           populations" (Numerical Simulation Study).
# Vignette: vignettes/fig02-bivariate-toy.qmd
#
# This is the single script for the bivariate toy: scenario builders
# (sourced by tests/testthat/helper.R) and the ADEMP pipeline.
# Set DECOVART_SOURCE_HELPERS=1 to load the functions without running
# the factorial grid.
#
# ── Factorial design ────────────────────────────────────────────────────────
#  Factor                  Levels
#  ─────────────────────── ───────────────────────────────────────────────────
#  Proportions             balanced (1/2, 1/2); highly unbalanced (19/20, 1/20)
#  Mean distance (CLD)     small: μ=(20,22)/(22,20); large: μ=(20,40)/(40,20)
#  Gene–gene corr CT1 (ρ)  −0.8 to +0.8, step 0.2  (9 levels)
#  Gene–gene corr CT2 (ρ)  −0.8 to +0.8, step 0.2  (9 levels)
#  Variance structure      homoscedastic (σ²=1,1); heteroscedastic (σ²=1,2)
#  Algorithms              NNLS, DeconRNASeq (LSEI), L-BFGS-B, gradient,
#                          Newton–Raphson, Marquardt–Levenberg, SA
#  ─────────────────────── ───────────────────────────────────────────────────
#  Total scenarios:  2 × 2 × 9 × 9 × 2 = 648
#  Replicates (n):   500 (full); 2 (smoke test)
#
# ── Usage ───────────────────────────────────────────────────────────────────
#  Full run:    Rscript scripts/fig02_bivariate_toy.R
#  Smoke test:  N_REPLICATES=2 Rscript scripts/fig02_bivariate_toy.R
#  Helpers:     DECOVART_SOURCE_HELPERS=1 source("scripts/fig02_bivariate_toy.R")
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/fig02/bivariate_benchmark.rds    – full benchmark list
#  output/fig02/bivariate_config.rds       – scenario config tibble
#  output/fig02/fig02_heatmap.pdf          – MSE heatmap (ρ₁ × ρ₂ grid)
#  output/fig02/fig02_raincloud.pdf        – raincloud of MC errors
#  output/fig02/fig02_forest.pdf           – forest / dot-whisker (ADEMP)
#  output/fig02/fig02_similarity.pdf       – algorithm similarity tile
###############################################################################

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
    "small unbalanced" = c(0.6, 0.4),
    "highly unbalanced" = c(0.05, 0.95)
  ),
  signature_matrices = list(
    "small OVL" = matrix(c(20, 40, 40, 20), nrow = 2)
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
      ID = paste0(
        "B",
        .data$scenario_idx,
        "_",
        ifelse(.data$variance == "homoscedastic", "Ho", "He")
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
        scenario_idx = scenario_idx,
        correlation_celltype1 = correlation_celltype1,
        correlation_celltype2 = correlation_celltype2,
        rho_ct1 = correlation_celltype1,
        rho_ct2 = correlation_celltype2,
        overlap = overlap,
        entropy = round(compute_shannon_entropy(p), digits = 3),
        proportions = proportion_name,
        proportion_name = proportion_name,
        variance = variance,
        centroids = centroids,
        centroid = centroids,
        true_theta = list(true_theta)
      )
    }
  ) |>
    dplyr::bind_rows()
}

#' Default deconvolution solvers for the bivariate toy benchmark
#'
#' Omits CIBERSORT (too few genes for nu-SVR tuning). Unconstrained OLS
#' (`lsfit`) is not shipped.
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


# Pipeline: skip when tests (or callers) only want the helpers.
if (!nzchar(Sys.getenv("DECOVART_SOURCE_HELPERS"))) {
  # ==========================================================================
  # SECTION 0 · Dependencies and paths
  # ==========================================================================

  if (!requireNamespace("DeCovarT", quietly = TRUE)) {
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::load_all(".", quiet = TRUE)
    } else {
      stop("Install DeCovarT or devtools before running this script.")
    }
  } else {
    library(DeCovarT)
  }

  stopifnot(
    requireNamespace("MixSim", quietly = TRUE),
    requireNamespace("ggplot2", quietly = TRUE)
  )

  OUT_DIR <- file.path("output", "fig02")
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

  N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "500"))
  if (!interactive() && N_REPL > 2L) {
    message(
      "fig02: running FULL benchmark (n = ",
      N_REPL,
      "). ",
      "Set N_REPLICATES=2 for a smoke test."
    )
  } else if (interactive()) {
    N_REPL <- 2L
    message("fig02: interactive session detected – using n = 2 (smoke test).")
  }

  SEED <- 20260903L
  set.seed(SEED)

  # ==========================================================================
  # SECTION 1 · GENERATIVE MODEL
  # ==========================================================================

  PROPORTIONS <- list(
    "balanced" = c(0.5, 0.5),
    "highly_unbalanced" = c(19 / 20, 1 / 20)
  )
  SIGNATURE_MATRICES <- list(
    "small_CLD" = matrix(
      c(20, 22, 22, 20),
      nrow = 2L,
      dimnames = list(NULL, NULL)
    ),
    "large_CLD" = matrix(
      c(20, 40, 40, 20),
      nrow = 2L,
      dimnames = list(NULL, NULL)
    )
  )
  CORR_SEQUENCE <- seq(-0.8, 0.8, by = 0.2)
  DIAGONAL_TERMS <- list(
    "homoscedastic" = c(1, 1),
    "heteroscedastic" = c(1, 2)
  )

  message("fig02 | building scenario config (648 scenarios) ...")
  scenario_config <- build_bivariate_scenario_config(
    proportions = PROPORTIONS,
    signature_matrices = SIGNATURE_MATRICES,
    corr_sequence = CORR_SEQUENCE,
    diagonal_terms = DIAGONAL_TERMS
  )
  message("fig02 | config built: ", nrow(scenario_config), " scenarios.")
  saveRDS(scenario_config, file.path(OUT_DIR, "bivariate_config.rds"))

  # ==========================================================================
  # SECTION 2 · INFERENCE
  # ==========================================================================

  ITMAX <- 200L
  EPSILON <- 1e-4
  deconvolution_functions <- bivariate_toy_deconvolution_functions(
    itmax = ITMAX,
    epsilon = EPSILON
  )

  message("fig02 | running benchmark (n = ", N_REPL, ") ...")
  bivariate_out <- run_simulation_benchmark(
    scenario_config = scenario_config,
    deconvolution_functions = deconvolution_functions,
    n = N_REPL,
    cores = 1L
  )
  message("fig02 | benchmark complete.")
  saveRDS(bivariate_out, file.path(OUT_DIR, "bivariate_benchmark.rds"))

  # ==========================================================================
  # SECTION 3 · VISUALISATIONS
  # ==========================================================================

  if (N_REPL < 10L) {
    message(
      "fig02 | skipping figure export (smoke-test run, n = ",
      N_REPL,
      ")."
    )
    message("fig02 | re-run with N_REPLICATES=500 to generate article figures.")
  } else {
    if (
      requireNamespace("ComplexHeatmap", quietly = TRUE) &&
        requireNamespace("circlize", quietly = TRUE) &&
        requireNamespace("viridis", quietly = TRUE)
    ) {
      heatmap_metrics <- dplyr::left_join(
        bivariate_out$optimisation,
        bivariate_out$config,
        by = intersect(
          names(bivariate_out$optimisation),
          names(bivariate_out$config)
        )
      )
      heatmap_list <- plot_correlation_Heatmap(
        distribution_metrics = heatmap_metrics,
        score_variable = "model_rmse"
      )
      pdf(file.path(OUT_DIR, "fig02_heatmap.pdf"), width = 10, height = 8)
      purrr::walk(heatmap_list, ComplexHeatmap::draw)
      grDevices::dev.off()
      message("fig02 | saved fig02_heatmap.pdf")
    } else {
      message("fig02 | ComplexHeatmap not available; skipping heatmap.")
    }

    if (requireNamespace("ggdist", quietly = TRUE)) {
      p_rain <- plot_mc_raincloud(
        bivariate_out,
        quantity = "error",
        facet_rows = "proportion_name",
        facet_cols = "centroid"
      )
      ggplot2::ggsave(
        file.path(OUT_DIR, "fig02_raincloud.pdf"),
        plot = p_rain,
        width = 12,
        height = 7
      )
      message("fig02 | saved fig02_raincloud.pdf")
    }

    p_forest <- plot_mc_forest(
      bivariate_out,
      facet_rows = "proportion_name",
      facet_cols = "centroid"
    )
    ggplot2::ggsave(
      file.path(OUT_DIR, "fig02_forest.pdf"),
      plot = p_forest,
      width = 12,
      height = 7
    )
    message("fig02 | saved fig02_forest.pdf")

    if (length(unique(bivariate_out$monte_carlo$algorithm)) >= 2L) {
      p_sim <- plot_algorithm_similarity(bivariate_out)
      ggplot2::ggsave(
        file.path(OUT_DIR, "fig02_similarity.pdf"),
        plot = p_sim,
        width = 6,
        height = 5
      )
      message("fig02 | saved fig02_similarity.pdf")
    }
  }

  message(
    "fig02 | done. Outputs in ",
    normalizePath(OUT_DIR, mustWork = FALSE)
  )
}
