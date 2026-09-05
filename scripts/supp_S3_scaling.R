###############################################################################
###############################################################################
###                                                                         ###
###     SUPPLEMENTARY S3 – HIGH-DIMENSIONAL SCALING GRID                   ###
###     J × G × cosine × condition × entropy · solver comparison            ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Background launch from the repository root (vanilla Rscript; no CLI
# parser — hyperparameters are hard-coded below). stdout and stderr go
# to logs/:
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore scripts/supp_S3_scaling.R \
#     > "logs/supp_S3_$(date +%F)_scaling.log" 2>&1 &
#
# Article:  DeCovarT – Supplementary: scalability
# Vignette: vignettes/fig03-variance-driven.qmd  {#sec-scenario-grid}
#           (pilot subset; full grid on HPC)
#
# ── Factorial design ─────────────────────────────────────────────────────────
#  Factor              Levels
#  ─────────────────── ─────────────────────────────────────────────────────
#  J (cell types)      3, 5, 10
#  G (genes)           50, 200, 500
#  Mean cosine         0.5, 0.8, 0.95
#  Cov. condition      low (u = 1.0); high (u = 0.1) [precision shift]
#  Proportions         balanced; highly unbalanced
#  Algorithms          NNLS, DeconRNASeq, Marquardt–Levenberg, L-BFGS-B
#  ─────────────────── ─────────────────────────────────────────────────────
#  Total:  3 × 3 × 3 × 2 × 2 = 108 scenario grids
#  Replicates (n):  100
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Rscript scripts/supp_S3_scaling.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/supp_S3/scaling_benchmark.rds
#  output/supp_S3/S3_metric_dots.pdf
###############################################################################

# ==============================================================================
# SECTION 0 · Dependencies and paths ----
# ==============================================================================

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

OUT_DIR <- file.path("output", "supp_S3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.ui_h1("Supplementary S3 · High-dimensional scaling")

N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "100"))

SEED <- 20260906L
set.seed(SEED)


# ==============================================================================
# SECTION 1 · GENERATIVE MODEL ----
#   Build scenario config for the (J, G, cosine, condition, entropy) grid.
# ==============================================================================

J_LEVELS <- c(3L, 5L, 10L)
G_LEVELS <- c(50L, 200L, 500L)
COSINE_LEVELS <- c(0.5, 0.8, 0.95)
PREC_SHIFTS <- c("low_condition" = 1.0, "high_condition" = 0.1)
PROP_CASES <- list(
  "balanced" = function(J) rep(1 / J, J),
  "highly_unbalanced" = function(J) {
    p <- rep(0.02, J)
    p[1L] <- 1 - 0.02 * (J - 1L)
    p
  }
)

full_grid <- tidyr::expand_grid(
  J = J_LEVELS,
  G = G_LEVELS,
  target_cosine = COSINE_LEVELS,
  condition_lbl = names(PREC_SHIFTS),
  prop_lbl = names(PROP_CASES)
)

.ui_info("Building {.val {nrow(full_grid)}} scenario configs.")

scenario_config_S3 <- purrr::pmap_dfr(
  full_grid,
  function(
    J,
    G,
    target_cosine,
    condition_lbl,
    prop_lbl
  ) {
    prec_shift <- PREC_SHIFTS[[condition_lbl]]
    p_fn <- PROP_CASES[[prop_lbl]]

    # Mean signature
    mu <- tryCatch(
      generate_mean_signature_matrix(
        n_genes = G,
        n_celltypes = J,
        mean_scale = 20,
        target_cosine = target_cosine,
        seed = SEED + G + J
      ),
      error = function(e) NULL
    )
    if (is.null(mu)) {
      return(NULL)
    }

    # Diagonal covariance (identity * shift, for scalability focus)
    Sigma <- array(0, dim = c(G, G, J))
    for (j in seq_len(J)) {
      diag(Sigma[,, j]) <- prec_shift
    }
    dimnames(Sigma) <- list(rownames(mu), rownames(mu), colnames(mu))

    p <- p_fn(J)

    tibble::tibble(
      J = J,
      G = G,
      target_cosine = target_cosine,
      condition_lbl = condition_lbl,
      prop_lbl = prop_lbl,
      entropy = round(compute_shannon_entropy(p), 3),
      true_theta = list(list(p = p, mu = mu, sigma = Sigma))
    )
  }
)

.ui_success(
  "Config: {.val {nrow(scenario_config_S3)}} valid scenarios."
)


# ==============================================================================
# SECTION 2 · INFERENCE ----
# ==============================================================================

deconvolution_functions_S3 <- list(
  "nnls" = list(FUN = deconvolute_ratios_nnls),
  "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
  "Marquardt-Levenberg" = list(
    FUN = deconvolute_ratios_Marquardt_Levenberg,
    additional_parameters = list(epsilon = 1e-4, itmax = 200L)
  ),
  "LBFGS" = list(
    FUN = deconvolute_ratios_L_BFGS_B,
    additional_parameters = list(epsilon = 1e-4, itmax = 200L)
  )
)

.ui_info("Running benchmark with {.val {N_REPL}} replicates.")
scaling_out <- run_simulation_benchmark(
  scenario_config = scenario_config_S3,
  deconvolution_functions = deconvolution_functions_S3,
  n = N_REPL,
  cores = 1L,
  verbose = TRUE
)
saveRDS(scaling_out, file.path(OUT_DIR, "scaling_benchmark.rds"))


# ==============================================================================
# SECTION 3 · VISUALISATIONS ----
# ==============================================================================

p_dots <- plot_mc_metric_dots(
  scaling_out,
  facet_rows = "condition_lbl",
  facet_cols = "target_cosine",
  metrics = c("rmse", "coverage")
)
ggplot2::ggsave(
  file.path(OUT_DIR, "S3_metric_dots.pdf"),
  plot = p_dots,
  width = 14,
  height = 8
)
.ui_success("Saved {.file S3_metric_dots.pdf}.")

.ui_success(
  "Done. Outputs in {.path {normalizePath(OUT_DIR, mustWork = FALSE)}}."
)
