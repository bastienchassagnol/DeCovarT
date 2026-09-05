###############################################################################
###############################################################################
###                                                                         ###
###     FIGURE 03 – VARIANCE-DRIVEN SCENARIO (article §2.2)                ###
###     J = 3 cell types · G = 50 genes · graph-constrained covariances     ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Background launch from the repository root (vanilla Rscript; no CLI
# parser — hyperparameters are hard-coded below). stdout and stderr go
# to logs/:
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore scripts/fig03_variance_driven.R \
#     > "logs/fig03_$(date +%F)_variance_driven.log" 2>&1 &
#
# Article:  DeCovarT – Section 2.2 (network topology separates mean-collinear
#           types). Vignette: vignettes/fig03-variance-driven.qmd
# Seed:     20260807
#
# ── Design ──────────────────────────────────────────────────────────────────
#  Means           one G x J signature; target Gram R with
#                  cos(μ1,μ2)=0.9, cos(μ1,μ3)=cos(μ2,μ3)=0.1
#  Graph models    SBM, Erdős–Rényi, hub/star  (independent per cell type)
#  Precision κ     well / moderate / ill  (precision_shift u)
#  Topology grid   3³ = 27 graph assignments
#  Covariance grid 27 x 3 = 81 (Σ_j, Ω_j) draws
#  Proportions     H*=1 (balanced); H*=0.5; H*=0.1
#  Algorithms      DeconRNASeq (LSEI), CIBERSORT, L-BFGS-B,
#                  Newton–Raphson, Marquardt–Levenberg
#                  (barycentre start; no NNLS)
#  Replicates (n)  50
#  Total scenarios 27 x 3 x 3 = 243
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Rscript scripts/fig03_variance_driven.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/fig03/hybrid_benchmark.rds
#  output/fig03/fig03_raincloud.pdf, fig03_forest.pdf, fig03_metric_dots.pdf
#  vignettes/figures/fig_network_topologies.png
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

if (!requireNamespace("igraph", quietly = TRUE)) {
  .ui_abort(
    "fig03 requires {.pkg igraph}. Install with {.code install.packages(\"igraph\")}."
  )
}
if (!requireNamespace("e1071", quietly = TRUE)) {
  .ui_abort(
    "fig03 requires {.pkg e1071} (CIBERSORT). Install with {.code install.packages(\"e1071\")}."
  )
}

OUT_DIR <- file.path("output", "fig03")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.ui_h1("Figure 03 · Variance-driven scenario")

N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "50"))
SEED <- 20260807L
set.seed(SEED)

N_GENES <- 50L
N_CELLTYPES <- 3L
MEAN_SCALE <- 10
ITMAX <- 200L
EPSILON <- 1e-4

celltype_palette <- c(
  celltype_1 = "#2B8C99",
  celltype_2 = "#F4B400",
  celltype_3 = "#C1443C"
)
graph_palette <- c(
  stochastic_block_model = "#4C72B0",
  erdos_renyi = "#999999",
  hub = "#C1443C"
)
ct_nms <- names(celltype_palette)

GRAPH_MODELS <- c("stochastic_block_model", "erdos_renyi", "hub")
# Spectral cushion u: lambda_min(Omega_j) = u after the shift, so
# kappa(Omega_j) = lambda_max / u. Smaller u -> worse-conditioned precision.
PRECISION_SHIFT <- c(well = 0.5, moderate = 0.1, ill = 0.02)
PRECISION_SCALE <- 0.3
PROP_INHIBITORY <- 0.5

GRAPH_PARAMS <- list(
  n_hubs = 1L,
  block_prob = c(0.4, 0.3, 0.3),
  p_within = 0.3,
  p_between = 0.01
)

TARGET_GRAM <- matrix(
  c(
    1,
    0.9,
    0.1,
    0.9,
    1,
    0.1,
    0.1,
    0.1,
    1
  ),
  nrow = 3L,
  dimnames = list(ct_nms, ct_nms)
)


# ==============================================================================
# SECTION 1 · GENERATIVE MODEL ----
#   Fixed mean signature (exact Gram). Factorial graphs x precision shift.
# ==============================================================================

mu <- generate_mean_signature_matrix(
  n_genes = N_GENES,
  n_celltypes = N_CELLTYPES,
  mean_scale = MEAN_SCALE,
  target_gram = TARGET_GRAM,
  seed = SEED,
  nonnegative = TRUE,
  celltype_names = ct_nms
)
cos_mu <- crossprod(mu)
norms <- sqrt(diag(cos_mu))
cos_mu <- cos_mu / tcrossprod(norms)
.ui_info(
  "Realised cosines: CT1–CT2 {.val {format(round(cos_mu[1, 2], 3), nsmall = 3)}}, CT1–CT3 {.val {format(round(cos_mu[1, 3], 3), nsmall = 3)}}, CT2–CT3 {.val {format(round(cos_mu[2, 3], 3), nsmall = 3)}}."
)

topology_grid <- tidyr::expand_grid(
  graph_ct1 = GRAPH_MODELS,
  graph_ct2 = GRAPH_MODELS,
  graph_ct3 = GRAPH_MODELS
)
kappa_grid <- tibble::tibble(
  kappa_label = names(PRECISION_SHIFT),
  precision_shift = unname(PRECISION_SHIFT)
)

cov_grid <- tidyr::expand_grid(topology_grid, kappa_grid)
.ui_info(
  "Covariance grid: {.val {nrow(topology_grid)}} topologies x {.val {nrow(kappa_grid)}} precision shifts = {.val {nrow(cov_grid)}}."
)

p_balanced <- composition_from_entropy(1, N_CELLTYPES, nms = ct_nms)
p_mod <- composition_from_entropy(0.5, N_CELLTYPES, nms = ct_nms)
p_rare <- composition_from_entropy(0.1, N_CELLTYPES, nms = ct_nms)
PROPORTIONS_3 <- list(
  "balanced" = p_balanced,
  "moderately unbalanced" = p_mod,
  "highly unbalanced" = p_rare
)
.ui_info(
  "Entropies: balanced {.val {round(compute_shannon_entropy(p_balanced), 3)}}, moderate {.val {round(compute_shannon_entropy(p_mod), 3)}}, high {.val {round(compute_shannon_entropy(p_rare), 3)}}."
)

.ui_info("Drawing one adjacency triple per topology (27 graphs).")
adj_by_topo <- purrr::pmap(
  topology_grid,
  function(graph_ct1, graph_ct2, graph_ct3) {
    models <- c(graph_ct1, graph_ct2, graph_ct3)
    lapply(seq_along(models), function(j) {
      withr::with_seed(
        SEED +
          j +
          10L * match(models[[j]], GRAPH_MODELS),
        generate_random_network_skeleton(
          n_genes = N_GENES,
          graph_model = models[[j]],
          graph_params = GRAPH_PARAMS
        )
      )
    })
  }
)

.ui_info("Completing signed precisions for each topology x kappa draw.")
cov_draws <- purrr::pmap(
  cov_grid,
  function(graph_ct1, graph_ct2, graph_ct3, kappa_label, precision_shift) {
    topo_idx <- which(
      topology_grid$graph_ct1 == graph_ct1 &
        topology_grid$graph_ct2 == graph_ct2 &
        topology_grid$graph_ct3 == graph_ct3
    )
    moments <- withr::with_seed(
      SEED +
        1000L * match(kappa_label, names(PRECISION_SHIFT)) +
        topo_idx,
      simulate_hierarchical_grn_moments(
        n_genes = N_GENES,
        n_celltypes = N_CELLTYPES,
        mean_scale = MEAN_SCALE,
        target_gram = TARGET_GRAM,
        precision_shift = precision_shift,
        precision_scale = PRECISION_SCALE,
        prop_inhibitory = PROP_INHIBITORY,
        nonnegative = TRUE,
        adjacency = adj_by_topo[[topo_idx]]
      )
    )
    list(
      sigma = moments$covariance_matrices,
      theta = moments$precision_matrices,
      adjacency = lapply(
        seq_len(N_CELLTYPES),
        function(j) moments$graph_structure$adjacency_matrices[,, j]
      )
    )
  }
)

scenario_config_3 <- purrr::map_dfr(
  seq_len(nrow(cov_grid)),
  function(i) {
    row <- cov_grid[i, , drop = FALSE]
    draw <- cov_draws[[i]]
    purrr::imap_dfr(PROPORTIONS_3, function(p, prop_name) {
      described <- describe_simulation_scenario(
        true_theta = list(
          p = p,
          mu = mu,
          sigma = draw$sigma,
          Theta = draw$theta
        ),
        adjacency = draw$adjacency
      )
      tibble::tibble(
        proportion_name = prop_name,
        entropy = round(compute_shannon_entropy(p), 3),
        graph_ct1 = row$graph_ct1,
        graph_ct2 = row$graph_ct2,
        graph_ct3 = row$graph_ct3,
        kappa_label = row$kappa_label,
        precision_shift = row$precision_shift,
        f_cov = described$descriptors$f_cov,
        kappa_sigma_p = described$descriptors$kappa_sigma_p,
        true_theta = list(list(
          p = p,
          mu = mu,
          sigma = draw$sigma,
          Theta = draw$theta,
          adjacency = draw$adjacency
        ))
      )
    })
  }
)
.ui_success(
  "Config built: {.val {nrow(scenario_config_3)}} scenarios."
)
saveRDS(scenario_config_3, file.path(OUT_DIR, "hybrid_config.rds"))


# ==============================================================================
# SECTION 2 · INFERENCE ----
# ==============================================================================

deconvolution_functions_3 <- list(
  "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
  "cibersort" = list(FUN = deconvolute_ratios_cibersort),
  "LBFGS" = list(
    FUN = deconvolute_ratios_L_BFGS_B,
    additional_parameters = list(
      epsilon = EPSILON,
      itmax = ITMAX,
      initial_p = "barycentre"
    )
  ),
  "Newton-Raphson" = list(
    FUN = deconvolute_ratios_Newton_Raphson,
    additional_parameters = list(
      epsilon = EPSILON,
      itmax = ITMAX,
      initial_p = "barycentre"
    )
  ),
  "Marquardt-Levenberg" = list(
    FUN = deconvolute_ratios_Marquardt_Levenberg,
    additional_parameters = list(
      epsilon = EPSILON,
      itmax = ITMAX,
      initial_p = "barycentre"
    )
  )
)

.ui_info(
  "Running ADEMP benchmark with {.val {N_REPL}} replicates."
)
hybrid_out <- run_simulation_benchmark(
  scenario_config = scenario_config_3,
  deconvolution_functions = deconvolution_functions_3,
  n = N_REPL,
  cores = 1L,
  verbose = TRUE
)
saveRDS(hybrid_out, file.path(OUT_DIR, "hybrid_benchmark.rds"))


# ==============================================================================
# SECTION 3 · VISUALISATIONS ----
# ==============================================================================

if (requireNamespace("ggdist", quietly = TRUE)) {
  p_rain <- plot_mc_raincloud(
    hybrid_out,
    quantity = "error",
    facet_rows = "proportion_name",
    facet_cols = "kappa_label"
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, "fig03_raincloud.pdf"),
    plot = p_rain,
    width = 12,
    height = 8
  )
  .ui_success("Saved {.file fig03_raincloud.pdf}.")
}

p_forest <- plot_mc_forest(
  hybrid_out,
  facet_rows = "proportion_name",
  facet_cols = "kappa_label"
)
ggplot2::ggsave(
  file.path(OUT_DIR, "fig03_forest.pdf"),
  plot = p_forest,
  width = 12,
  height = 8
)
.ui_success("Saved {.file fig03_forest.pdf}.")

p_dots <- plot_mc_metric_dots(
  hybrid_out,
  facet_rows = "proportion_name",
  facet_cols = "kappa_label",
  metrics = c("rmse", "mae", "coverage")
)
ggplot2::ggsave(
  file.path(OUT_DIR, "fig03_metric_dots.pdf"),
  plot = p_dots,
  width = 12,
  height = 8
)
.ui_success("Saved {.file fig03_metric_dots.pdf}.")

# One panel per graph generator (same G, independent of the factorial).
static_network_path <- file.path(
  "vignettes",
  "figures",
  "fig_network_topologies.png"
)
dir.create(
  dirname(static_network_path),
  recursive = TRUE,
  showWarnings = FALSE
)
grDevices::png(static_network_path, width = 2700, height = 900, res = 220)
graphics::par(mfrow = c(1L, 3L), mar = c(1, 1, 3, 1))
set.seed(SEED)
for (gm in GRAPH_MODELS) {
  adj <- generate_random_network_skeleton(
    n_genes = N_GENES,
    graph_model = gm,
    graph_params = GRAPH_PARAMS
  )
  graph <- igraph::graph_from_adjacency_matrix(
    adj,
    mode = "undirected",
    diag = FALSE
  )
  igraph::V(graph)$color <- graph_palette[[gm]]
  igraph::plot.igraph(
    graph,
    vertex.size = 4,
    vertex.label = NA,
    vertex.frame.color = "#2f3e4f",
    edge.color = "#4a5560",
    edge.width = 0.8,
    layout = igraph::layout_with_fr(graph),
    main = gm
  )
}
grDevices::dev.off()
.ui_success("Wrote {.path {static_network_path}}.")

.ui_success(
  "Done. Outputs in {.path {normalizePath(OUT_DIR, mustWork = FALSE)}}."
)
