###############################################################################
###############################################################################
###                                                                         ###
###     SUPPLEMENTARY S1 – IDENTIFIABILITY EDGE CASES                      ###
###     Boundary zeros · Same-mean identifiability · CI method comparison   ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Article:  DeCovarT – Appendix S1 (regular-case MLE + identifiability)
# Vignette: vignettes/supp-S1-identifiability.qmd
# Theory:   vignettes/theory-DeCovarT-MLE-properties.qmd
#
# ── Design ──────────────────────────────────────────────────────────────────
#  Regular interior (live in the vignette; repeated here as a smoke check)
#    affine equivariance; bulk perturbation; Dirichlet ILR starts
#  S1a: boundary null (p₃=0)  J=3, G=10        Wald / profile / bootstrap CI
#  S1b: same means / diff Σ   J=2, G=5         ρ₁₂ ∈ {0, 0.5, 0.9}
#  S1c: multimodal lik.        J=2, G=2         ρ ∈ {−0.8, 0, 0.8}; low CLD
#
#  Replicates (n): 1000 (full); 2 (smoke test)
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Full run:    Rscript scripts/supp_S1_identifiability.R
#  Smoke test:  N_REPLICATES=2 Rscript scripts/supp_S1_identifiability.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/supp_S1/S1a_boundary.rds
#  output/supp_S1/S1b_same_means.rds
#  output/supp_S1/S1c_multimodal.rds
###############################################################################

# ==============================================================================
# SECTION 0 · Dependencies and paths
# ==============================================================================

if (!requireNamespace("DeCovarT", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("Install DeCovarT or devtools before running this script.")
  }
} else {
  library(DeCovarT)
}

OUT_DIR <- file.path("output", "supp_S1")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "1000"))
if (interactive()) {
  N_REPL <- 2L
  message("supp_S1: interactive session – smoke test (n = 2).")
}

SEED <- 20260904L
set.seed(SEED)


# ==============================================================================
# SECTION 1 · GENERATIVE MODEL
#   Three sub-cases designed to stress identifiability limits.
# ==============================================================================

# ── S1a: Active boundary null – p₃ = 0 ──────────────────────────────────────
#   J = 3, G = 10.  True composition: (0.6, 0.4, 0.0).
#   Test whether Wald, profile-likelihood, and parametric-bootstrap CIs
#   achieve nominal coverage near a simplex face.
J_1a <- 3L
G_1a <- 10L
set.seed(SEED + 1L)
mu_1a <- generate_mean_signature_matrix(
  n_genes = G_1a,
  n_celltypes = J_1a,
  mean_scale = 20,
  target_cosine = 0.3,
  seed = SEED + 1L
)
Sigma_1a <- array(
  0,
  dim = c(G_1a, G_1a, J_1a),
  dimnames = list(rownames(mu_1a), rownames(mu_1a), colnames(mu_1a))
)
for (j in seq_len(J_1a)) {
  diag(Sigma_1a[,, j]) <- 1
}
p_1a <- c(0.6, 0.4, 0.0) # zero boundary

config_1a <- tibble::tibble(
  sub_case = "S1a_boundary_null",
  p3_true = p_1a[3L],
  true_theta = list(list(p = p_1a, mu = mu_1a, sigma = Sigma_1a))
)
message("supp_S1 | S1a generative model built (J=3, G=10, p₃=0).")


# ── S1b: Same means, different Σ – identifiability from covariance only ──────
#   J = 2, G = 5.  μ₁ = μ₂ (identical columns).  Σ₁ and Σ₂ differ.
#   Vary the off-diagonal correlation ρ₁₂ so that the covariance signal grows.
G_1b <- 5L
rho_grid_1b <- c(0, 0.5, 0.9)
mu_1b <- matrix(
  20,
  nrow = G_1b,
  ncol = 2L,
  dimnames = list(paste0("g", seq_len(G_1b)), c("ct1", "ct2"))
)
p_1b <- c(0.5, 0.5)

config_1b <- purrr::map_dfr(rho_grid_1b, function(rho) {
  R <- matrix(rho, G_1b, G_1b)
  diag(R) <- 1
  Sigma_b <- array(
    0,
    dim = c(G_1b, G_1b, 2L),
    dimnames = list(rownames(mu_1b), rownames(mu_1b), colnames(mu_1b))
  )
  Sigma_b[,, 1L] <- R # CT1: equicorrelation ρ
  diag(Sigma_b[,, 2L]) <- 1 # CT2: diagonal (independent)
  tibble::tibble(
    sub_case = "S1b_same_means",
    rho_ct1 = rho,
    true_theta = list(list(p = p_1b, mu = mu_1b, sigma = Sigma_b))
  )
})
message(
  "supp_S1 | S1b generative model built (",
  nrow(config_1b),
  " ρ levels)."
)


# ── S1c: Multimodal likelihood – low CLD, symmetric composition ──────────────
#   J = 2, G = 2.  Nearly identical means → near-flat likelihood ridge.
#   Correlation ρ ∈ {−0.8, 0, 0.8} alters ridge shape.
mu_1c <- matrix(
  c(20, 21, 21, 20),
  nrow = 2L,
  dimnames = list(c("g1", "g2"), c("ct1", "ct2"))
)
rho_grid_1c <- c(-0.8, 0, 0.8)
p_1c <- c(0.5, 0.5)

config_1c <- purrr::map_dfr(rho_grid_1c, function(rho) {
  R <- matrix(rho, 2L, 2L)
  diag(R) <- 1
  Sigma_c <- array(
    c(R, R),
    dim = c(2L, 2L, 2L),
    dimnames = list(c("g1", "g2"), c("g1", "g2"), c("ct1", "ct2"))
  )
  tibble::tibble(
    sub_case = "S1c_multimodal",
    rho = rho,
    true_theta = list(list(p = p_1c, mu = mu_1c, sigma = Sigma_c))
  )
})
message(
  "supp_S1 | S1c generative model built (",
  nrow(config_1c),
  " ρ levels)."
)


# ── Regular interior smoke check (matches vignette live chunks) ─────────────
#   Cheap: one bulk, three Dirichlet starts. Always runs, even when
#   N_REPL is a smoke-test value.
genes_int <- paste0("g", 1:3)
cts_int <- paste0("ct", 1:2)
mu_int <- matrix(
  c(20, 40, 25, 35, 30, 22),
  nrow = 3,
  dimnames = list(genes_int, cts_int)
)
Sigma_int <- array(
  0,
  dim = c(3, 3, 2),
  dimnames = list(genes_int, genes_int, cts_int)
)
Sigma_int[,, 1L] <- diag(c(1.0, 1.2, 0.8))
Sigma_int[,, 2L] <- diag(c(1.5, 0.9, 1.1))
p_int <- c(0.6, 0.4)
set.seed(SEED)
y_int <- drop(mu_int %*% p_int + stats::rnorm(3, sd = 0.2))
starts_int <- lapply(c(11L, 22L, 33L), function(s) {
  set.seed(s)
  starting_simplex(2L, "dirichlet", nms = cts_int)
})
p_ml_starts <- vapply(
  starts_int,
  function(p0) {
    deconvolute_ratios_Marquardt_Levenberg(
      y_int,
      mu_int,
      Sigma_int,
      itmax = 80L,
      initial_p = p0
    )
  },
  numeric(2)
)
start_spread <- max(apply(p_ml_starts, 1L, stats::sd))
message(
  "supp_S1 | interior ILR-start spread (max SD across starts) = ",
  signif(start_spread, 3),
  if (start_spread < 0.02) " (solvers agree)." else " (WARNING: spread)."
)
saveRDS(
  list(starts = starts_int, marquardt = p_ml_starts, spread = start_spread),
  file.path(OUT_DIR, "S1_interior_starts.rds")
)


# ==============================================================================
# SECTION 2 · INFERENCE
#   S1a: compare Wald vs profile vs bootstrap CI coverage near boundary.
#   S1b / S1c: Marquardt–Levenberg only (focus is on likelihood shape).
# ==============================================================================

deconv_fns_simple <- list(
  "Marquardt-Levenberg" = list(
    FUN = deconvolute_ratios_Marquardt_Levenberg,
    additional_parameters = list(epsilon = 1e-4, itmax = 200L)
  )
)

message("supp_S1 | running S1a (n = ", N_REPL, ") ...")
out_1a <- run_simulation_benchmark(
  scenario_config = config_1a,
  deconvolution_functions = deconv_fns_simple,
  n = N_REPL,
  cores = 1L
)
saveRDS(out_1a, file.path(OUT_DIR, "S1a_boundary.rds"))

message("supp_S1 | running S1b (n = ", N_REPL, ") ...")
out_1b <- run_simulation_benchmark(
  scenario_config = config_1b,
  deconvolution_functions = deconv_fns_simple,
  n = N_REPL,
  cores = 1L
)
saveRDS(out_1b, file.path(OUT_DIR, "S1b_same_means.rds"))

message("supp_S1 | running S1c (n = ", N_REPL, ") ...")
out_1c <- run_simulation_benchmark(
  scenario_config = config_1c,
  deconvolution_functions = deconv_fns_simple,
  n = N_REPL,
  cores = 1L
)
saveRDS(out_1c, file.path(OUT_DIR, "S1c_multimodal.rds"))
message("supp_S1 | inference complete.")


# ==============================================================================
# SECTION 3 · VISUALISATIONS
#   S1a: forest plot of CI widths (Wald vs profile); coverage table.
#   S1b: log-likelihood contour vs ρ facet.
#   S1c: log-likelihood 1D slice per ρ.
# ==============================================================================

if (N_REPL < 10L) {
  message("supp_S1 | skipping figure export (smoke-test run).")
} else {
  # S1b: dot-metric comparison across ρ levels
  p_dots_1b <- plot_mc_metric_dots(
    out_1b,
    facet_cols = "rho_ct1",
    metrics = c("rmse", "bias", "coverage")
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, "S1b_metric_dots.pdf"),
    plot = p_dots_1b,
    width = 10,
    height = 5
  )

  # S1c: raincloud of errors per ρ
  if (requireNamespace("ggdist", quietly = TRUE)) {
    p_rain_1c <- plot_mc_raincloud(
      out_1c,
      quantity = "error",
      facet_cols = "rho"
    )
    ggplot2::ggsave(
      file.path(OUT_DIR, "S1c_raincloud.pdf"),
      plot = p_rain_1c,
      width = 10,
      height = 5
    )
  }
  message("supp_S1 | figures saved.")
}

message("supp_S1 | done. Outputs in ", normalizePath(OUT_DIR, mustWork = FALSE))
