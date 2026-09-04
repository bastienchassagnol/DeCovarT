###############################################################################
###############################################################################
###                                                                         ###
###     SUPPLEMENTARY S4 – IMPACT OF COVARIANCE MODELLING                  ###
###     OLS vs diagonal vs global GLS vs full Σ_j (DeCovarT)               ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Article:  DeCovarT – Supplementary: covariance modelling comparison
# Vignette: vignettes/supp-S4-covariance-modeling.qmd
#           {#sec-perturb} extends the two-covariance comparison to a
#           fuller three-way grid.
#
# ── Factorial design ─────────────────────────────────────────────────────────
#  Factor              Levels
#  ─────────────────── ─────────────────────────────────────────────────────
#  True Σ structure    isotropic (σ²I_G); block-diagonal; full (Markov net)
#  Inference model     (a) NNLS / ignore Σ
#                      (b) CT-diagonal Σ (zero off-diagonals per CT)
#                      (c) global GLS (Σ_global = (1/J²) Σ_j, known at p=1/J)
#                      (d) full Σ_j (DeCovarT, Marquardt–Levenberg)
#  J                   3 (fixed)
#  G                   50 (fixed; re-uses fig03 network seed)
#  Mean cosine         0.7, 0.95
#  Replicates (n)      200 (full); 2 (smoke test)
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Smoke test:  N_REPLICATES=2 Rscript scripts/supp_S4_covariance_modeling.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/supp_S4/covariance_modeling_benchmark.rds
#  output/supp_S4/S4_forest.pdf
#  output/supp_S4/S4_metric_dots.pdf
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

OUT_DIR <- file.path("output", "supp_S4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "200"))
if (interactive()) {
  N_REPL <- 2L
  message("supp_S4: interactive – smoke test (n = 2).")
}

SEED <- 20260907L
set.seed(SEED)

J_FIXED <- 3L
G_FIXED <- 50L


# ==============================================================================
# SECTION 1 · GENERATIVE MODEL
#   Three true covariance structures crossed with two mean cosine levels.
# ==============================================================================

COSINE_LEVELS <- c(0.7, 0.95)
p_eq <- rep(1 / J_FIXED, J_FIXED)

# Helper: build a sparse block-diagonal Sigma from igraph-based precision
.make_grn_sigma <- function(G, J, cosine, seed) {
  mu <- generate_mean_signature_matrix(
    n_genes       = G,
    n_celltypes   = J,
    mean_scale    = 20,
    target_cosine = cosine,
    seed          = seed
  )
  Omega_arr <- array(0, dim = c(G, G, J),
                     dimnames = list(rownames(mu), rownames(mu), colnames(mu)))
  for (j in seq_len(J)) {
    set.seed(seed + j)
    adj   <- generate_random_network_skeleton(G, "erdos_renyi",
                                             list(p_erdos_renyi = 0.2))
    w_adj <- assign_iid_signed_weights(adj, prop_inhibitory = 0.5,
                                       weight_magnitude = 0.4)
    Omega_arr[, , j] <- build_normalised_precision(w_adj, precision_shift = 1)
  }
  list(mu = mu, Sigma = build_covariance_array_from_precision(Omega_arr))
}

message("supp_S4 | building generative models ...")
scenario_config_S4 <- purrr::map_dfr(COSINE_LEVELS, function(rho) {
  moments_iso <- {
    mu_r <- generate_mean_signature_matrix(G_FIXED, J_FIXED, 20, rho,
                                           seed = SEED + round(rho * 100))
    Sigma_iso <- array(0, dim = c(G_FIXED, G_FIXED, J_FIXED),
                       dimnames = list(rownames(mu_r), rownames(mu_r), colnames(mu_r)))
    for (j in seq_len(J_FIXED)) diag(Sigma_iso[, , j]) <- 1
    list(mu = mu_r, Sigma = Sigma_iso)
  }

  moments_grn <- .make_grn_sigma(G_FIXED, J_FIXED, rho, SEED + 50L)

  purrr::map_dfr(
    list(
      "isotropic" = moments_iso,
      "GRN_based" = moments_grn
    ),
    function(m) {
      tibble::tibble(
        target_cosine = rho,
        true_sigma_structure = c("isotropic", "GRN_based")[
          match(list(m), list(moments_iso, moments_grn))
        ],
        true_theta = list(list(
          p     = p_eq,
          mu    = m$mu,
          sigma = m$Sigma
        ))
      )
    },
    .id = "true_sigma_structure"
  )
})
# Fix the .id column (purrr sets it from map names)
scenario_config_S4 <- scenario_config_S4 |>
  dplyr::select(-dplyr::any_of("true_sigma_structure...3")) |>
  dplyr::distinct()

message("supp_S4 | config: ", nrow(scenario_config_S4), " scenarios.")


# ==============================================================================
# SECTION 2 · INFERENCE
#   Four solvers: NNLS, CT-diagonal ML, global GLS, full DeCovarT.
#   The GLS solver uses fixed_gls_covariance() at the equal composition.
# ==============================================================================

# CT-diagonal solver: zero off-diagonals of each Σ_j before passing to ML
.deconv_ct_diagonal <- function(y, mu, Sigma, ...) {
  Sigma_diag <- Sigma
  for (j in seq_len(dim(Sigma)[3L])) {
    Sigma_diag[, , j] <- diag(diag(Sigma[, , j]))
  }
  deconvolute_ratios_Marquardt_Levenberg(y, mu, Sigma_diag, ...)
}

# Global GLS solver: fixed covariance at equi-balanced composition
.deconv_global_gls <- function(y, mu, Sigma, ...) {
  p_equal <- rep(1 / dim(Sigma)[3L], dim(Sigma)[3L])
  w_global <- fixed_gls_covariance(Sigma, p = p_equal, diagonal = FALSE)
  deconvolute_ratios_gls(y, mu, w_global)
}

deconvolution_functions_S4 <- list(
  "nnls" = list(FUN = deconvolute_ratios_nnls),
  "ct_diagonal_ML" = list(
    FUN = .deconv_ct_diagonal,
    additional_parameters = list(epsilon = 1e-4, itmax = 200L)
  ),
  "global_gls" = list(FUN = .deconv_global_gls),
  "full_DeCovarT" = list(
    FUN = deconvolute_ratios_Marquardt_Levenberg,
    additional_parameters = list(epsilon = 1e-4, itmax = 200L)
  )
)

message("supp_S4 | running benchmark (n = ", N_REPL, ") ...")
cov_modeling_out <- run_simulation_benchmark(
  scenario_config        = scenario_config_S4,
  deconvolution_functions = deconvolution_functions_S4,
  n                      = N_REPL,
  cores                  = 1L)
saveRDS(cov_modeling_out, file.path(OUT_DIR, "covariance_modeling_benchmark.rds"))
message("supp_S4 | benchmark complete.")


# ==============================================================================
# SECTION 3 · VISUALISATIONS
# ==============================================================================

if (N_REPL < 10L) {
  message("supp_S4 | skipping figures (smoke-test run).")
} else {
  p_forest <- plot_mc_forest(
    cov_modeling_out,
    facet_rows = "true_sigma_structure",
    facet_cols = "target_cosine"
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, "S4_forest.pdf"),
    plot   = p_forest,
    width  = 12,
    height = 7
  )

  p_dots <- plot_mc_metric_dots(
    cov_modeling_out,
    facet_rows = "true_sigma_structure",
    facet_cols = "target_cosine",
    metrics    = c("rmse", "bias", "coverage")
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, "S4_metric_dots.pdf"),
    plot   = p_dots,
    width  = 14,
    height = 7
  )
  message("supp_S4 | figures saved.")
}

message("supp_S4 | done. Outputs in ", normalizePath(OUT_DIR, mustWork = FALSE))
