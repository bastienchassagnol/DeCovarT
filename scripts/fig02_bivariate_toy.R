###############################################################################
###                                                                         ###
###           FIGURE 02 - BIVARIATE TOY MODEL (article Section 2.2)         ###
###       J = 2 cell types, G = 2 genes, full factorial design              ###
###                                                                         ###
###############################################################################
#
# Run from the repository root; stdout and stderr go to logs/:
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore scripts/fig02_bivariate_toy.R \
#     > "logs/fig02_$(date +%F)_bivariate_toy.log" 2>&1 &
#
# Run size (environment variables; see
# scripts/benchmark_metrics/benchmark_artefacts.R):
#   DECOVART_N_REPLICATES  Monte Carlo bulks per scenario   (default 200)
#   DECOVART_N_BOOT        bootstrap refits per bulk          (default 100)
#   DECOVART_CORES         workers for the per-sample loop    (default 1)
#                          (workers load the installed DeCovarT, so run
#                          devtools::install() before a parallel run)
#   DECOVART_POSTPROCESS_ONLY=1 redraws the figures from the saved benchmark.
# A dry run: DECOVART_N_REPLICATES=5 DECOVART_N_BOOT=5 Rscript ...
#
# Factorial design
#  Proportions             balanced (1/2, 1/2); moderately unbalanced
#                          (17/20, 3/20); highly unbalanced (99/100, 1/100)
#  Centroid gap            small: mu = (20, 22) / (22, 20);
#                          large: mu = (20, 40) / (40, 20)
#  Correlation, type 1     -0.8 to 0.8 by 0.2 (9 levels)
#  Correlation, type 2     -0.8 to 0.8 by 0.2 (9 levels)
#  Variances               homoscedastic (1, 1); heteroscedastic (1, 2)
#  Scenarios               3 x 2 x 9 x 9 x 2 = 972
#  Scenario ID             B{index}_{Ho|He}_{Ba|Mo|Hi}_{Sm|Lg}
#
# Methods
#  lsei            DeconRNASeq's simplex-constrained least squares, also the
#                  point estimate of DECALS (no public code)
#  NNLS            non-negative least squares with a free scale, renormalised
#                  afterwards; it differs from lsei, which keeps the sum-to-one
#                  equality and the bounds at every iteration
#  Newton-Raphson  DeCovarT, trust-region Newton in ILR coordinates
#  L-BFGS-B        DeCovarT, L-BFGS-B on the closure of the orthant
#  Firth           DeCovarT with the Jeffreys penalty
#  CIBERSORT (whose nu-SVR would regress two bulk values on two signatures
#  and an intercept) and the second-generation methods (MEAD, RNA-Sieve,
#  ReDeconv, whose variance estimates need many more genes than cell types)
#  are run in fig03 and fig05 only.
#
# Confidence intervals at every DeCovarT estimate: Wald (delta method), ILR
# (Wald in ILR coordinates mapped back), profile likelihood, parametric
# bootstrap refitted by the same optimiser, and the oracle half-width at the
# true composition (Newton-Raphson only). Firth: Wald and ILR at the
# penalised estimate.
#
# Outputs (output/fig02/)
#  fig02_benchmark.rds   settings, design, truth, every metric table
#  fig02_plot_data.rds   the data behind every figure, keyed by file stem
#  density_visualisations/*.pdf, performance_visualisations/*.pdf
#
# Sourcing this file from a test or a vignette defines the builders only:
# the pipeline runs for Rscript scripts/fig02_bivariate_toy.R or an
# interactive session.
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
    proportions = names(proportion_list),
    correlation_celltype1 = corr_sequence,
    correlation_celltype2 = corr_sequence,
    variance = names(diagonal_terms)
  ) |>
    dplyr::mutate(
      scenario_idx = dplyr::row_number(),
      ID = .encode_bivariate_id(
        .data$scenario_idx,
        .data$variance,
        .data$proportions,
        .data$centroids
      )
    )

  purrr::pmap(
    design,
    function(
      centroids,
      proportions,
      correlation_celltype1,
      correlation_celltype2,
      variance,
      scenario_idx,
      ID
    ) {
      mu <- signature_matrices[[centroids]]
      p <- proportion_list[[proportions]]
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
        proportions = proportions,
        variance = variance,
        centroids = centroids,
        true_theta = list(true_theta)
      )
    }
  ) |>
    dplyr::bind_rows()
}

#' Methods and confidence intervals of the bivariate benchmark
#'
#' @param itmax,epsilon Iteration cap and convergence tolerance of the
#'   DeCovarT solvers.
#' @param n_boot Parametric-bootstrap refits per bulk; `NULL` drops the
#'   intervals (tests).
#'
#' @return Named list for [run_simulation_benchmark()].
bivariate_toy_deconvolution_functions <- function(
  itmax = 10L,
  epsilon = 1e-3,
  n_boot = NULL
) {
  solver <- list(epsilon = epsilon, itmax = itmax, initial_p = "barycentre")
  interval_menu <- function(methods, optimiser = "newton") {
    if (is.null(n_boot)) {
      return(NULL)
    }
    list(methods = methods, n_boot = n_boot, optimiser = optimiser)
  }
  list(
    "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
    "NNLS" = list(FUN = deconvolute_ratios_nnls),
    "Newton-Raphson" = list(
      FUN = deconvolute_ratios_Newton_Raphson,
      additional_parameters = solver,
      ci = interval_menu(c("wald", "ilr", "profile", "bootstrap", "oracle"))
    ),
    "L-BFGS-B" = list(
      FUN = deconvolute_ratios_L_BFGS_B,
      additional_parameters = solver,
      ci = interval_menu(c("wald", "ilr", "profile", "bootstrap"), "lbfgsb")
    ),
    "Firth" = list(
      FUN = fit_firth_decovart,
      ci = interval_menu(c("wald", "ilr"))
    )
  )
}

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
  devtools::load_all(".", quiet = TRUE)
  source(file.path("scripts", "benchmark_metrics", "benchmark_artefacts.R"))
  source(file.path("scripts", "visualisation", "load_visualisation.R"))
  stopifnot(requireNamespace("MixSim", quietly = TRUE))

  OUT_DIR <- file.path("output", "fig02")
  BENCHMARK_FILE <- file.path(OUT_DIR, "fig02_benchmark.rds")
  SEED <- 20260903L
  ITMAX <- 200L
  EPSILON <- 1e-4
  CORES <- as.integer(Sys.getenv("DECOVART_CORES", "1"))
  settings <- benchmark_settings(default_replicates = 200L, default_boot = 100L)

  .ui_h1("Figure 02: bivariate toy model")

  if (!settings$postprocess_only) {
    set.seed(SEED)
    scenario_config <- build_bivariate_scenario_config()
    .ui_info(
      "{.val {nrow(scenario_config)}} scenarios, N = {.val {settings$n_replicates}}, B = {.val {settings$n_boot}}."
    )
    started <- proc.time()[["elapsed"]]
    bivariate_out <- run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = bivariate_toy_deconvolution_functions(
        itmax = ITMAX,
        epsilon = EPSILON,
        n_boot = settings$n_boot
      ),
      n = settings$n_replicates,
      cores = CORES,
      verbose = TRUE,
      return_data = TRUE
    )
    .ui_success(
      "Benchmark finished in {.val {round(proc.time()[['elapsed']] - started)}} s."
    )
    save_benchmark(
      bivariate_out,
      scenario_config,
      benchmark_meta(SEED, settings),
      BENCHMARK_FILE
    )
  }

  benchmark <- read_benchmark(BENCHMARK_FILE)
  plot_data <- figure_plot_data(benchmark, figure = "fig02")
  saveRDS(plot_data, file.path(OUT_DIR, "fig02_plot_data.rds"))
  render_figure_books(plot_data, OUT_DIR)
  .ui_success("Figure 02 written to {.file {OUT_DIR}}.")
}
