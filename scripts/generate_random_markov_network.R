# ---- generate_random_markov_network.R --------------------------------------
##
## End-to-end synthetic benchmark verifying the DeCovarT feature-selection
## pipeline on a deliberately hard 3-cell-type / 50-gene scenario:
##
##   * cell types 1 and 2 are strongly correlated at the MEAN level (they
##     cannot be told apart from mu alone) and are instead discriminated
##     through their PRECISION MATRIX / network topology (a BLGGM-style
##     hybrid design: hub/stochastic-block vs star/key-driver);
##   * cell type 3 is set apart by a block of highly specific marker genes
##     wired as a scale-free network;
##   * a block of genes is equally expressed (and equally wired) across all
##     three types, i.e. carries no discriminative signal at all.
##
## The script has six numbered sections, each independently re-runnable
## after the previous ones:
##   1) mean signature matrix  mu (50 genes x 3 cell types)
##   2) cell-type-specific block-structured topologies and precisions
##   3) N = 200 simulated bulk mixtures (+ latent cell-type-specific tensor)
##   4) feature-selection pipeline (pre-selection, then NSGA-II refinement)
##   5) gene-level summary table across pipeline stages
##   6) visuals: ComplexHeatmap of mu, interactive visNetwork topologies
##
## Run from the package root: Rscript scripts/generate_random_markov_network.R

stop_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing package '",
      pkg,
      "'. Install with e.g. pak::pak(\"",
      pkg,
      "\").",
      call. = FALSE
    )
  }
}
invisible(lapply(
  c("mco", "visNetwork", "ComplexHeatmap", "circlize", "readr"),
  stop_if_missing
))

devtools::load_all(".", quiet = TRUE)
set.seed(20260807)

out_dir <- file.path(getwd(), "output", "random_markov_network")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# A single colour chart reused throughout (mean-profile figures, heatmap,
# network plots) so that cell-type identity stays visually consistent.
celltype_palette <- c(
  celltype_1 = "#2B8C99", # teal   - hub / SBM topology
  celltype_2 = "#F4B400", # gold   - star / key-driver topology
  celltype_3 = "#C1443C" # coral  - scale-free marker topology
)
block_palette <- c(
  shared_12_vs_3 = "#4C72B0", # DE(12 vs 3), but "DM" between 1 and 2
  marker_3 = "#C1443C", # DE marker of cell type 3
  equal_all = "#999999" # EE: equally expressed / equally wired
)

# ==============================================================================
# 1. Mean signature matrix mu (G = 50 genes x J = 3 cell types)
# ==============================================================================
## Design (50 genes total):
##   * 30 "shared_12_vs_3" genes: near-identical mean in cell types 1 & 2
##     (high target cosine), background level in cell type 3. These genes
##     cannot separate type 1 from type 2 at the mean level -> the network
##     topology built in Section 2 does that job instead.
##   * 10 "marker_3" genes: near-orthogonal marker block, highly expressed
##     only in cell type 3, poorly expressed in types 1 and 2.
##   * 10 "equal_all" genes: identical mean in all three cell types (no
##     discriminative signal whatsoever; "EE" in the scDD / muscat sense).
##
## generate_mean_signature_matrix() is called twice, iteratively:
##   (i)  on the 30 "shared_12_vs_3" genes, for the pair (cell type 1,
##        cell type 2), with a HIGH target cosine (strongly correlated);
##   (ii) on a temporary 20-gene pool, for the pair (merged {1,2}, cell
##        type 3), with a LOW target cosine (near-orthogonal); only the
##        10 genes whose *private* block belongs to cell type 3 are kept,
##        so that both flanking blocks (private-to-merged and
##        private-to-type-3) are respected and cell types 1 and 2 receive
##        the same (near-zero) background level on this marker block.
## The 10 "equal_all" genes need no cosine machinery (their columns must be
## identical by construction) and are assigned a flat baseline directly.

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
  rep(c("shared_12_vs_3", "marker_3", "equal_all"), c(n_g12, n_g3, n_geq)),
  levels = c("shared_12_vs_3", "marker_3", "equal_all")
)
names(gene_block) <- gene_names

mean_scale <- 10
target_cosine_12 <- 0.95 # (i) cell types 1 & 2: strongly correlated means
target_cosine_3 <- 0.05 # (ii) merged {1,2} vs type 3: near orthogonal
background_level <- 0.5 # low baseline for "off-block" genes
equal_level <- 2 # flat baseline for the equally-expressed block, of the
# same order of magnitude as the shared_12_vs_3 block so that the heatmap
# colour scale reflects "no discriminative signal" rather than "highest
# expression level"

## --- (i) 30 genes: cell type 1 vs cell type 2, high target cosine ----------
mu_12 <- generate_mean_signature_matrix(
  n_genes = n_g12,
  n_celltypes = 2L,
  mean_scale = mean_scale,
  target_cosine = target_cosine_12,
  gene_names = gene_names[idx_g12],
  celltype_names = c("celltype_1", "celltype_2")
)

## --- (ii) 10 genes: merged {1,2} vs cell type 3, near-orthogonal -----------
## Ask for twice as many genes as needed (2 * n_g3) so that both the
## "private to merged {1,2}" and "private to cell type 3" blocks have
## exactly n_g3 genes each; keep only the latter block.
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

## --- Assemble the full G x J mean signature matrix -------------------------
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

mean_objectives <- compute_mean_profile_objectives(mu)
message(
  "Section 1 - mean profile objectives: mean |cosine| = ",
  signif(mean_objectives$mean_abs_cosine, 3),
  ", sum Euclidean distance = ",
  signif(mean_objectives$sum_euclidean_distance, 3)
)

# ==============================================================================
# 2. Cell-type-specific block-structured topologies (BLGGM-style hybrid)
# ==============================================================================
## Each cell type gets its OWN 50 x 50 precision matrix via
## simulate_hierarchical_grn_moments(..., adjacency = ...). Following the
## BLGGM hybrid design (global gene blocks + local topology per block),
## each cell type's adjacency is block-diagonal: the *local* random-graph
## model differs by (cell type, gene block) pair, with no edges between
## blocks:
##
##   gene block        | cell type 1        | cell type 2   | cell type 3
##   -------------------|---------------------|----------------|-------------
##   shared_12_vs_3 (30)| stochastic_block_model (hub-like, cascading  | hub, n_hubs = 1 (star, single | erdos_renyi (background;
##                      | pathway modules)                             | key-driver gene)              | uninformative here)
##   marker_3       (10)| erdos_renyi (background)                    | erdos_renyi (background)      | scale_free (cell type 3's
##                      |                                              |                                | expressed marker genes)
##   equal_all      (10)| erdos_renyi                                 | erdos_renyi                    | erdos_renyi
##
## i.e. cell types 1 and 2 differ only on the mean-collinear block (30
## genes), which is exactly where the mean signature could not separate
## them; the "equal_all" block is wired identically (ER) in every cell
## type, matching its complete lack of discriminative signal.

precision_shift <- 0.1
precision_scale <- 0.3
prop_inhibitory <- 0.5 # "as many positive [inhibitory] as activatory edges"

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
      graph_params = list(block_prob = c(0.4, 0.3, 0.3), p_within = 0.3)
    ),
    list(idx = idx_g3, graph_model = "erdos_renyi", graph_params = list()),
    list(idx = idx_geq, graph_model = "erdos_renyi", graph_params = list())
  ),
  celltype_2 = list(
    list(
      idx = idx_g12,
      graph_model = "hub",
      graph_params = list(n_hubs = 1L) # single star: one key-driver gene
    ),
    list(idx = idx_g3, graph_model = "erdos_renyi", graph_params = list()),
    list(idx = idx_geq, graph_model = "erdos_renyi", graph_params = list())
  ),
  celltype_3 = list(
    list(idx = idx_g12, graph_model = "erdos_renyi", graph_params = list()),
    list(
      idx = idx_g3,
      graph_model = "scale_free",
      graph_params = list(power = 1, edges_per_node = 1L)
    ),
    list(idx = idx_geq, graph_model = "erdos_renyi", graph_params = list())
  )
)

adjacency_list <- lapply(
  celltype_names,
  function(ct) build_block_adjacency(n_genes, block_defs_by_celltype[[ct]])
)
names(adjacency_list) <- celltype_names

## Package generator: per-cell-type Omega_j / Sigma_j from the hybrid
## adjacencies. Mean profiles from Section 1 overwrite the default cosine
## construction (which cannot encode the EE / DE / DM block design).
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
adjacency_matrices <- lapply(seq_len(3L), function(j) {
  adjacency_list[[j]]
})
names(adjacency_matrices) <- celltype_names
weighted_adjacencies <- lapply(seq_len(3L), function(j) {
  w <- grn_moments$graph_structure$weighted_adjacencies[,, j]
  dimnames(w) <- list(gene_names, gene_names)
  w
})
names(weighted_adjacencies) <- celltype_names

message(
  "Section 2 - built ",
  length(adjacency_matrices),
  " cell-type-specific ",
  n_genes,
  " x ",
  n_genes,
  " precision matrices (edges per type: ",
  paste(
    vapply(adjacency_matrices, function(a) sum(a) / 2L, numeric(1)),
    collapse = ", "
  ),
  ")."
)

# ==============================================================================
# 3. Simulate N = 200 bulk mixtures (+ latent cell-type-specific tensor)
# ==============================================================================
## simulate_bulk_mixture() draws, for every sample i = 1, ..., N, latent
## purified profiles x_{.j}^(i) ~ N(mu_{.j}, Sigma_j) independently for each
## cell type j, then forms the bulk sample as the p-weighted convolution.
## It already returns both the bulk matrix Y and the full latent tensor
## X in $latent_profiles (G x J x N); no extra plumbing is required.

n_samples <- 200L
p_true <- c(celltype_1 = 0.5, celltype_2 = 0.3, celltype_3 = 0.2)

bulk_simulation <- simulate_bulk_mixture(
  signature_matrix = mu,
  Sigma = Sigma,
  p = p_true,
  n = n_samples
)
latent_profiles <- bulk_simulation$latent_profiles # G x J x N
bulk_Y <- bulk_simulation$Y # G x N

message(
  "Section 3 - simulated ",
  n_samples,
  " bulk profiles; latent tensor dim = ",
  paste(dim(latent_profiles), collapse = " x "),
  "."
)

# ==============================================================================
# 4. Feature-selection pipeline
# ==============================================================================
true_theta <- list(p = p_true, mu = mu, sigma = Sigma)
check_true_theta(true_theta)

## ---- 4.i Pre-selection of the gene universe -------------------------------
## Two complementary, purely univariate/marginal screens:
##  a) compute_average_overlap(), applied gene-by-gene to the TRUE moments
##     (marginal mean mu[g, ] and marginal variance Sigma[g, g, ] under a
##     3-component Gaussian mixture with weights p_true): a *low* average
##     overlap flags a well-separated gene.
##  b) compute_glmnet_gene_scores(), applied once to the full latent tensor
##     from Section 3: a *large* multinomial elastic-net coefficient flags
##     a gene with strong multivariate class signal.
## The gene universe G0 is the union of the top `keep_fraction` genes under
## each criterion (a permissive pre-screen; the Pareto stage below tightens
## the panel further).

per_gene_overlap <- vapply(
  seq_len(n_genes),
  function(g) {
    theta_g <- list(
      p = p_true,
      mu = matrix(mu[g, ], nrow = 1L),
      sigma = Sigma[g, g, , drop = FALSE]
    )
    compute_average_overlap(theta_g)
  },
  numeric(1)
)
names(per_gene_overlap) <- gene_names

glmnet_scores <- compute_glmnet_gene_scores(
  expression_profiles = latent_profiles,
  celltype_labels = celltype_names,
  alpha = 0.5
)

keep_fraction <- 0.7
n_keep <- ceiling(keep_fraction * n_genes)
keep_by_overlap <- gene_names[order(per_gene_overlap)][seq_len(n_keep)]
keep_by_glmnet <- gene_names[order(-glmnet_scores)][seq_len(n_keep)]
gene_universe <- union(keep_by_overlap, keep_by_glmnet)

message(
  "Section 4.i - pre-selected gene universe G0: ",
  length(gene_universe),
  " / ",
  n_genes,
  " genes (",
  sum(gene_block[gene_universe] == "equal_all"),
  " of which are 'equal_all')."
)

## ---- 4.ii Multi-objective (NSGA-II) refinement of G0 ----------------------
## Bi-objective panel search over binary gene masks on G0 (@sec-pareto-panel,
## vignettes/feature-selection.qmd), minimising:
##  a) the panel condition number kappa_2(mu_{G}) (stability);
##  b) MINUS the average pairwise Jeffreys divergence between the 3 cell
##     types (compute_average_jeffreys()). NSGA-II always MINIMISES; since a
##     *larger* Jeffreys divergence means better-separated cell types, the
##     fitness returns "-Jeffreys" so that minimising it maximises actual
##     class separation, as recommended in the vignette ("maximise average
##     Jeffreys, if preferred").
n_universe <- length(gene_universe)
penalty <- 1e6

panel_objectives <- function(x) {
  active <- gene_universe[x > 0.5]
  if (length(active) < 3L) {
    return(c(penalty, penalty))
  }
  mu_active <- mu[active, , drop = FALSE]
  sigma_active <- Sigma[active, active, , drop = FALSE]
  kappa_val <- tryCatch(
    kappa(mu_active, exact = TRUE),
    error = function(e) penalty
  )
  jeffreys_val <- tryCatch(
    compute_average_jeffreys(list(
      mu = mu_active,
      sigma = sigma_active,
      p = p_true
    )),
    error = function(e) 0
  )
  c(kappa_val, -jeffreys_val)
}

nsga2_result <- mco::nsga2(
  fn = panel_objectives,
  idim = n_universe,
  odim = 2L,
  lower.bounds = rep(0, n_universe),
  upper.bounds = rep(1, n_universe),
  popsize = 80L,
  generations = 60L
)

## Best-compromise solution on the final Pareto front: the non-dominated,
## non-degenerate point closest (min-max normalised) to the ideal point.
pareto_idx <- which(
  nsga2_result$pareto.optimal &
    nsga2_result$value[, 1L] < penalty &
    nsga2_result$value[, 2L] < penalty
)
pareto_values <- nsga2_result$value[pareto_idx, , drop = FALSE]
pareto_pars <- nsga2_result$par[pareto_idx, , drop = FALSE]
value_range <- pmax(apply(pareto_values, 2L, function(v) diff(range(v))), 1e-8)
normalised_values <- sweep(
  sweep(pareto_values, 2L, apply(pareto_values, 2L, min), `-`),
  2L,
  value_range,
  `/`
)
best_compromise_idx <- which.min(sqrt(rowSums(normalised_values^2)))
best_mask <- pareto_pars[best_compromise_idx, ] > 0.5
final_panel <- gene_universe[best_mask]

message(
  "Section 4.ii - Pareto front: ",
  length(pareto_idx),
  " non-dominated panels; best compromise panel size = ",
  length(final_panel),
  " (kappa = ",
  signif(pareto_values[best_compromise_idx, 1L], 3),
  ", avg. Jeffreys = ",
  signif(-pareto_values[best_compromise_idx, 2L], 3),
  ")."
)

# ==============================================================================
# 5. Gene-level summary table across pipeline stages
# ==============================================================================
## Connects gene blocks to the scDD (Korthauer et al. 2016) / muscat
## (Crowell et al. 2020) differential-distribution taxonomy: "equal_all" is
## an EE (equivalent expression) block; "marker_3" is a DE block (cell type
## 3 vs {1, 2}); "shared_12_vs_3" behaves as DE against cell type 3 but as
## a DM (differential modality) block *between* cell types 1 and 2, since
## the mean is shared and only the network topology (hub/SBM vs star)
## differs (see vignettes/figures/fig_dd_taxonomy_ee_de_dm.png).
dd_category <- c(
  shared_12_vs_3 = "DE (vs type 3) / DM (type 1 vs type 2)",
  marker_3 = "DE (type 3 marker)",
  equal_all = "EE (equivalent expression)"
)

gene_selection_summary <- tibble::tibble(
  gene = gene_names,
  block = as.character(gene_block),
  dd_category = dd_category[as.character(gene_block)],
  overlap_score = per_gene_overlap,
  glmnet_score = glmnet_scores,
  in_universe_g0 = gene_names %in% gene_universe,
  in_final_panel = gene_names %in% final_panel
)

readr::write_csv(
  gene_selection_summary,
  file.path(out_dir, "gene_selection_summary.csv")
)
true_moments_rds <- list(
  p = p_true,
  mu = mu,
  X = latent_profiles,
  Sigma = Sigma,
  Theta = Theta
)
saveRDS(
  true_moments_rds,
  file.path(out_dir, "true_grn_moments.rds")
)
message(
  "Section 5 - wrote gene-level summary table (",
  nrow(gene_selection_summary),
  " genes) and true-moment RDS (p, mu, X, Sigma, Theta): ",
  out_dir
)
print(
  gene_selection_summary |>
    dplyr::count(.data$block, .data$in_universe_g0, .data$in_final_panel)
)

# ==============================================================================
# 6. Visuals
# ==============================================================================

## ---- 6.i ComplexHeatmap of mu, rows split by gene block -------------------
heatmap_path <- file.path(out_dir, "heatmap_mean_profiles.png")
grDevices::png(heatmap_path, width = 1400, height = 1800, res = 200)
mean_heatmap <- ComplexHeatmap::Heatmap(
  mu,
  name = "mean\nexpression",
  col = circlize::colorRamp2(
    c(0, mean(mu), max(mu)),
    c("#F7FBFF", "#6BAED6", "#08306B")
  ),
  row_split = gene_block,
  cluster_row_slices = FALSE,
  row_title_rot = 0,
  cluster_columns = FALSE,
  column_labels = celltype_names,
  top_annotation = ComplexHeatmap::HeatmapAnnotation(
    `cell type` = celltype_names,
    col = list(`cell type` = celltype_palette),
    show_legend = TRUE
  ),
  left_annotation = ComplexHeatmap::rowAnnotation(
    `gene block` = gene_block,
    col = list(`gene block` = block_palette),
    show_legend = TRUE
  ),
  show_row_names = FALSE,
  row_title_gp = grid::gpar(fontsize = 10),
  heatmap_legend_param = list(title = "mean\nexpression")
)
ComplexHeatmap::draw(mean_heatmap)
grDevices::dev.off()
message("Section 6.i - wrote: ", heatmap_path)

## ---- 6.ii Interactive visNetwork topology, one widget per cell type -------
build_visnetwork <- function(adjacency, weighted_adjacency, celltype) {
  degree <- rowSums(adjacency)
  nodes <- data.frame(
    id = gene_names,
    label = gene_names,
    group = as.character(gene_block),
    value = pmax(degree, 1),
    color.background = unname(block_palette[as.character(gene_block)]),
    color.border = "#2f3e4f",
    title = paste0(
      "<b>",
      gene_names,
      "</b><br/>block: ",
      as.character(gene_block),
      "<br/>degree: ",
      degree
    ),
    stringsAsFactors = FALSE
  )

  edge_idx <- which(upper.tri(adjacency) & adjacency != 0L, arr.ind = TRUE)
  edge_weight <- weighted_adjacency[edge_idx]
  edges <- data.frame(
    from = gene_names[edge_idx[, 1L]],
    to = gene_names[edge_idx[, 2L]],
    color = ifelse(edge_weight > 0, "#C1443C", "#2B8C99"),
    title = ifelse(edge_weight > 0, "inhibitory", "activatory"),
    stringsAsFactors = FALSE
  )

  legend_nodes <- data.frame(
    label = paste0("block: ", names(block_palette)),
    shape = "dot",
    color = unname(block_palette),
    stringsAsFactors = FALSE
  )

  visNetwork::visNetwork(
    nodes = nodes,
    edges = edges,
    main = paste0(celltype, " precision-matrix topology"),
    width = "100%",
    height = "800px"
  ) |>
    visNetwork::visLegend(
      useGroups = FALSE,
      addNodes = legend_nodes,
      addEdges = data.frame(
        label = c("inhibitory (+)", "activatory (-)"),
        color = c("#C1443C", "#2B8C99"),
        stringsAsFactors = FALSE
      )
    ) |>
    visNetwork::visInteraction(navigationButtons = TRUE) |>
    visNetwork::visOptions(highlightNearest = TRUE) |>
    visNetwork::visPhysics(stabilization = TRUE)
}

for (ct in celltype_names) {
  network_plot <- build_visnetwork(
    adjacency_matrices[[ct]],
    weighted_adjacencies[[ct]],
    ct
  )
  html_path <- file.path(out_dir, paste0("network_", ct, ".html"))
  visNetwork::visSave(network_plot, html_path)
  message("Section 6.ii - wrote: ", html_path)
}

## ---- 6.iii Static network figure, for inclusion in the pkgdown vignette ---
## The interactive visNetwork widgets above (Section 6.ii) are exploratory
## artefacts (large htmlwidgets-dependency folders) and stay untracked under
## output/. This single static PNG is the minimal, git-tracked counterpart
## that vignettes/synthetic-scenarios.qmd embeds and pkgdown renders.
static_network_path <- file.path(
  "vignettes",
  "figures",
  "fig_network_topologies.png"
)
grDevices::png(static_network_path, width = 2700, height = 1150, res = 220)
graphics::layout(
  matrix(c(1L, 2L, 3L, 4L, 4L, 4L), nrow = 2L, byrow = TRUE),
  heights = c(5, 1)
)
graphics::par(mar = c(1, 1, 3, 1))
set.seed(20260807)
for (ct in celltype_names) {
  graph <- igraph::graph_from_adjacency_matrix(
    adjacency_matrices[[ct]],
    mode = "undirected",
    diag = FALSE
  )
  edge_ends <- igraph::ends(graph, igraph::E(graph), names = FALSE)
  edge_weight <- weighted_adjacencies[[ct]][
    cbind(edge_ends[, 1L], edge_ends[, 2L])
  ]
  igraph::V(graph)$color <- unname(block_palette[as.character(gene_block)])
  igraph::plot.igraph(
    graph,
    vertex.size = 4,
    vertex.label = NA,
    vertex.frame.color = "#2f3e4f",
    edge.color = ifelse(edge_weight > 0, "#C1443C", "#2B8C99"),
    edge.width = 1.2,
    layout = igraph::layout_with_fr(graph),
    main = paste0(
      ct,
      " (",
      sub("celltype_", "type ", ct),
      ")"
    )
  )
}
graphics::par(mar = c(0, 0, 0, 0))
graphics::plot.new()
graphics::legend(
  "center",
  legend = c(names(block_palette), "inhibitory (+)", "activatory (-)"),
  col = c(unname(block_palette), "#C1443C", "#2B8C99"),
  pch = c(19, 19, 19, NA, NA),
  lty = c(NA, NA, NA, 1, 1),
  lwd = 2,
  bty = "n",
  horiz = TRUE
)
grDevices::dev.off()
message("Section 6.iii - wrote: ", static_network_path)
