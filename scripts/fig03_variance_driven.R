###############################################################################
###############################################################################
###                                                                         ###
###     FIGURE 03 – VARIANCE-DRIVEN HYBRID SCENARIO (article §2.2)         ###
###     J = 3 cell types · G = 50 genes · GRN-based covariances             ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Article:  DeCovarT – Section 2.2 (network topology separates mean-collinear
#           types). Vignette: vignettes/fig03-variance-driven.qmd
# Seed:     20260807 (canonical; matches vignette geometry tables)
#
# Single script for this figure: mean signature, block-structured
# precisions, ADEMP benchmark, and the static network PNG.
# Feature-selection (NSGA-II) is Appendix S6, not this scenario.
#
# ── Design ──────────────────────────────────────────────────────────────────
#  Gene block        Genes  CT 1 topology  CT 2 topology  CT 3 topology
#  ─────────────     ─────  ─────────────  ─────────────  ─────────────
#  shared_12_vs_3    30     SBM / hub      star           Erdős–Rényi
#  marker_3          10     Erdős–Rényi    Erdős–Rényi    scale-free
#  equal_all         10     Erdős–Rényi    Erdős–Rényi    Erdős–Rényi
#
#  Proportions     balanced (1/3,1/3,1/3); mod. unbalanced (0.5,0.3,0.2);
#                  highly unbalanced (0.7,0.2,0.1)
#  Algorithms      NNLS, DeconRNASeq (LSEI), Marquardt–Levenberg
#  Replicates (n)  200 (full); 2 (smoke test)
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Full run:    Rscript scripts/fig03_variance_driven.R
#  Smoke test:  N_REPLICATES=2 Rscript scripts/fig03_variance_driven.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  data/synthetic_networks/true_grn_moments.rds
#  output/fig03/hybrid_benchmark.rds
#  output/fig03/fig03_raincloud.pdf, fig03_forest.pdf, fig03_metric_dots.pdf
#  vignettes/figures/fig_network_topologies.png
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

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop(
    "fig03 requires igraph. Install with install.packages(\"igraph\")."
  )
}

OUT_DIR <- file.path("output", "fig03")
DATA_DIR <- file.path("data", "synthetic_networks")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

N_REPL <- as.integer(Sys.getenv("N_REPLICATES", "200"))
if (interactive()) {
  N_REPL <- 2L
  message("fig03: interactive session – smoke test (n = 2).")
}

SEED <- 20260807L
set.seed(SEED)

celltype_palette <- c(
  celltype_1 = "#2B8C99",
  celltype_2 = "#F4B400",
  celltype_3 = "#C1443C"
)
block_palette <- c(
  shared_12_vs_3 = "#4C72B0",
  marker_3 = "#C1443C",
  equal_all = "#999999"
)


# ==============================================================================
# SECTION 1 · GENERATIVE MODEL
#   Canonical G = 50, J = 3 hybrid: types 1–2 mean-collinear on 30 genes;
#   type 3 has a marker block; 10 genes are a null (equal_all) block.
# ==============================================================================

GRN_CACHE <- file.path(DATA_DIR, "true_grn_moments.rds")

if (file.exists(GRN_CACHE)) {
  message("fig03 | loading cached GRN moments from ", GRN_CACHE)
  true_grn_moments <- readRDS(GRN_CACHE)
} else {
  message("fig03 | generating GRN moments (seed ", SEED, ") ...")

  n_g12 <- 30L
  n_g3 <- 10L
  n_geq <- 10L
  n_genes <- n_g12 + n_g3 + n_geq
  celltype_names <- names(celltype_palette)

  gene_names <- c(
    paste0("gene_12_", seq_len(n_g12)),
    paste0("gene_3_", seq_len(n_g3)),
    paste0("gene_eq_", seq_len(n_geq))
  )
  idx_g12 <- seq_len(n_g12)
  idx_g3 <- n_g12 + seq_len(n_g3)
  idx_geq <- n_g12 + n_g3 + seq_len(n_geq)

  gene_block <- factor(
    rep(
      c("shared_12_vs_3", "marker_3", "equal_all"),
      c(n_g12, n_g3, n_geq)
    ),
    levels = c("shared_12_vs_3", "marker_3", "equal_all")
  )
  names(gene_block) <- gene_names

  mean_scale <- 10
  target_cosine_12 <- 0.95
  target_cosine_3 <- 0.05
  background_level <- 0.5
  equal_level <- 2

  mu_12 <- generate_mean_signature_matrix(
    n_genes = n_g12,
    n_celltypes = 2L,
    mean_scale = mean_scale,
    target_cosine = target_cosine_12,
    gene_names = gene_names[idx_g12],
    celltype_names = c("celltype_1", "celltype_2")
  )
  mu_3_pool <- generate_mean_signature_matrix(
    n_genes = 2L * n_g3,
    n_celltypes = 2L,
    mean_scale = mean_scale,
    target_cosine = target_cosine_3,
    gene_names = paste0("pool_", seq_len(2L * n_g3)),
    celltype_names = c("celltype_12_merged", "celltype_3")
  )
  mu_3_block <- mu_3_pool[(n_g3 + 1L):(2L * n_g3), , drop = FALSE]
  rownames(mu_3_block) <- gene_names[idx_g3]

  mu <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = 3L,
    dimnames = list(gene_names, celltype_names)
  )
  mu[idx_g12, "celltype_1"] <- mu_12[, "celltype_1"]
  mu[idx_g12, "celltype_2"] <- mu_12[, "celltype_2"]
  mu[idx_g12, "celltype_3"] <- background_level
  mu[idx_g3, "celltype_1"] <- mu_3_block[, "celltype_12_merged"]
  mu[idx_g3, "celltype_2"] <- mu_3_block[, "celltype_12_merged"]
  mu[idx_g3, "celltype_3"] <- mu_3_block[, "celltype_3"]
  mu[idx_geq, ] <- equal_level
  if (any(mu < 0)) {
    mu <- mu + abs(min(mu))
  }

  precision_shift <- 0.1
  precision_scale <- 0.3
  prop_inhibitory <- 0.5

  build_block_adjacency <- function(n_genes, block_defs) {
    adjacency <- matrix(0L, n_genes, n_genes)
    for (blk in block_defs) {
      adjacency[blk$idx, blk$idx] <- generate_random_network_skeleton(
        n_genes = length(blk$idx),
        graph_model = blk$graph_model,
        graph_params = blk$graph_params
      )
    }
    dimnames(adjacency) <- list(gene_names, gene_names)
    adjacency
  }

  block_defs_by_celltype <- list(
    celltype_1 = list(
      list(
        idx = idx_g12,
        graph_model = "stochastic_block_model",
        graph_params = list(
          block_prob = c(0.4, 0.3, 0.3),
          p_within = 0.3
        )
      ),
      list(
        idx = idx_g3,
        graph_model = "erdos_renyi",
        graph_params = list()
      ),
      list(
        idx = idx_geq,
        graph_model = "erdos_renyi",
        graph_params = list()
      )
    ),
    celltype_2 = list(
      list(
        idx = idx_g12,
        graph_model = "hub",
        graph_params = list(n_hubs = 1L)
      ),
      list(
        idx = idx_g3,
        graph_model = "erdos_renyi",
        graph_params = list()
      ),
      list(
        idx = idx_geq,
        graph_model = "erdos_renyi",
        graph_params = list()
      )
    ),
    celltype_3 = list(
      list(
        idx = idx_g12,
        graph_model = "erdos_renyi",
        graph_params = list()
      ),
      list(
        idx = idx_g3,
        graph_model = "scale_free",
        graph_params = list(power = 1, edges_per_node = 1L)
      ),
      list(
        idx = idx_geq,
        graph_model = "erdos_renyi",
        graph_params = list()
      )
    )
  )

  adjacency_list <- lapply(
    celltype_names,
    function(ct) {
      build_block_adjacency(n_genes, block_defs_by_celltype[[ct]])
    }
  )
  names(adjacency_list) <- celltype_names

  grn_moments <- simulate_hierarchical_grn_moments(
    n_genes = n_genes,
    n_celltypes = 3L,
    mean_scale = mean_scale,
    target_cosine = 0,
    precision_shift = precision_shift,
    precision_scale = precision_scale,
    prop_inhibitory = prop_inhibitory,
    adjacency = adjacency_list
  )
  grn_moments$mean_profiles <- mu
  rownames(grn_moments$mean_profiles) <- gene_names
  colnames(grn_moments$mean_profiles) <- celltype_names

  Sigma <- grn_moments$covariance_matrices
  Theta <- grn_moments$precision_matrices
  dimnames(Sigma) <- list(gene_names, gene_names, celltype_names)
  dimnames(Theta) <- list(gene_names, gene_names, celltype_names)

  true_grn_moments <- list(
    p = c(celltype_1 = 0.5, celltype_2 = 0.3, celltype_3 = 0.2),
    mu = mu,
    Sigma = Sigma,
    Theta = Theta,
    gene_block = gene_block,
    adjacency_list = adjacency_list,
    weighted_adjacencies = lapply(seq_len(3L), function(j) {
      w <- grn_moments$graph_structure$weighted_adjacencies[,, j]
      dimnames(w) <- list(gene_names, gene_names)
      w
    }),
    seed = SEED
  )
  names(true_grn_moments$weighted_adjacencies) <- celltype_names
  saveRDS(true_grn_moments, GRN_CACHE)
  message("fig03 | GRN moments saved to ", GRN_CACHE)
}

mu <- true_grn_moments$mu
Sigma <- true_grn_moments$Sigma
ct_nms <- colnames(mu)
gene_block <- true_grn_moments$gene_block
adjacency_list <- true_grn_moments$adjacency_list
weighted_adjacencies <- true_grn_moments$weighted_adjacencies

cos12 <- as.numeric(
  crossprod(mu[, 1L], mu[, 2L]) /
    (sqrt(sum(mu[, 1L]^2)) * sqrt(sum(mu[, 2L]^2)))
)
message(sprintf("fig03 | pairwise cosine CT1–CT2: %.3f", cos12))


# ==============================================================================
# SECTION 2 · INFERENCE
# ==============================================================================

PROPORTIONS_3 <- list(
  "balanced" = c(1 / 3, 1 / 3, 1 / 3),
  "mod_unbalanced" = c(0.5, 0.3, 0.2),
  "highly_unbalanced" = c(0.7, 0.2, 0.1)
)

scenario_config_3 <- purrr::imap_dfr(PROPORTIONS_3, function(p, prop_name) {
  tibble::tibble(
    proportion_name = prop_name,
    entropy = round(compute_shannon_entropy(p), 3),
    true_theta = list(list(
      p = stats::setNames(p, ct_nms),
      mu = mu,
      sigma = Sigma
    ))
  )
})

ITMAX <- 200L
EPSILON <- 1e-4
deconvolution_functions_3 <- list(
  "nnls" = list(FUN = deconvolute_ratios_nnls),
  "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
  "Marquardt-Levenberg" = list(
    FUN = deconvolute_ratios_Marquardt_Levenberg,
    additional_parameters = list(epsilon = EPSILON, itmax = ITMAX)
  )
)

message("fig03 | running ADEMP benchmark (n = ", N_REPL, ") ...")
hybrid_out <- run_simulation_benchmark(
  scenario_config = scenario_config_3,
  deconvolution_functions = deconvolution_functions_3,
  n = N_REPL,
  cores = 1L
)
message("fig03 | benchmark complete.")
saveRDS(hybrid_out, file.path(OUT_DIR, "hybrid_benchmark.rds"))


# ==============================================================================
# SECTION 3 · VISUALISATIONS
# ==============================================================================

if (N_REPL < 10L) {
  message("fig03 | skipping ADEMP figure export (smoke-test run).")
} else {
  if (requireNamespace("ggdist", quietly = TRUE)) {
    p_rain <- plot_mc_raincloud(
      hybrid_out,
      quantity = "error",
      facet_rows = "proportion_name"
    )
    ggplot2::ggsave(
      file.path(OUT_DIR, "fig03_raincloud.pdf"),
      plot = p_rain,
      width = 12,
      height = 7
    )
    message("fig03 | saved fig03_raincloud.pdf")
  }

  p_forest <- plot_mc_forest(
    hybrid_out,
    facet_rows = "proportion_name"
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, "fig03_forest.pdf"),
    plot = p_forest,
    width = 10,
    height = 6
  )
  message("fig03 | saved fig03_forest.pdf")

  p_dots <- plot_mc_metric_dots(
    hybrid_out,
    facet_rows = "proportion_name",
    metrics = c("rmse", "mae", "coverage")
  )
  ggplot2::ggsave(
    file.path(OUT_DIR, "fig03_metric_dots.pdf"),
    plot = p_dots,
    width = 12,
    height = 6
  )
  message("fig03 | saved fig03_metric_dots.pdf")
}

# Static network PNG for the vignette (always, even on smoke tests).
if (!is.null(adjacency_list) && !is.null(gene_block)) {
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
  grDevices::png(static_network_path, width = 2700, height = 1150, res = 220)
  graphics::layout(
    matrix(c(1L, 2L, 3L, 4L, 4L, 4L), nrow = 2L, byrow = TRUE),
    heights = c(5, 1)
  )
  graphics::par(mar = c(1, 1, 3, 1))
  set.seed(SEED)
  gene_names <- rownames(mu)
  for (ct in ct_nms) {
    graph <- igraph::graph_from_adjacency_matrix(
      adjacency_list[[ct]],
      mode = "undirected",
      diag = FALSE
    )
    edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
    edge_weight <- weighted_adjacencies[[ct]][
      cbind(edge_ends[, 1L], edge_ends[, 2L])
    ]
    igraph::V(graph)$color <- unname(
      block_palette[as.character(gene_block)]
    )
    igraph::plot.igraph(
      graph,
      vertex.size = 4,
      vertex.label = NA,
      vertex.frame.color = "#2f3e4f",
      edge.color = ifelse(edge_weight > 0, "#C1443C", "#2B8C99"),
      edge.width = 1.2,
      layout = igraph::layout_with_fr(graph),
      main = paste0(ct, " (", sub("celltype_", "type ", ct), ")")
    )
  }
  graphics::par(mar = c(0, 0, 0, 0))
  graphics::plot.new()
  graphics::legend(
    "center",
    legend = c(
      names(block_palette),
      "inhibitory (+)",
      "activatory (-)"
    ),
    col = c(unname(block_palette), "#C1443C", "#2B8C99"),
    pch = c(19, 19, 19, NA, NA),
    lty = c(NA, NA, NA, 1, 1),
    lwd = 2,
    bty = "n",
    horiz = TRUE
  )
  grDevices::dev.off()
  message("fig03 | wrote ", static_network_path)
}

message("fig03 | done. Outputs in ", normalizePath(OUT_DIR, mustWork = FALSE))
