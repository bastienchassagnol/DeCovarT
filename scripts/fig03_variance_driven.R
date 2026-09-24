###############################################################################
###############################################################################
###                                                                         ###
###     FIGURE 03 – COVARIANCE-DRIVEN SCENARIO (article §2.2)              ###
###     J = 3 cell types · G = 20 genes · graph-constrained covariances     ###
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
#     > "logs/fig03_$(date +%F)_covariance_driven.log" 2>&1 &
#
# Redraw figures from finished hybrid_*.rds (does not refit ADEMP):
#   FIG03_POSTPROCESS_ONLY=1 Rscript --no-save --no-restore \
#     scripts/fig03_variance_driven.R
#
# Article:  DeCovarT – Section 2.2 (network topology separates mean-collinear
#           types). Vignette: vignettes/fig03-covariance-driven.qmd
# Seed:     20260807
#
# ── Design ──────────────────────────────────────────────────────────────────
#  Means           one G x J signature; target Gram R with
#                  cos(μ1,μ2)=0.9, cos(μ1,μ3)=cos(μ2,μ3)=0.1
#  Graph models    scale-free (Barabasi-Albert) and cluster /
#                  stochastic block model.
#  Topology grid   2^2 = 4 assignments on mean-collinear CT1 and CT2;
#                  CT3 (mean-separated) is held at scale-free
#  Overlap         low / moderate / high target MixSim BarOmega
#                  (scale Sigma_j, keep precision zeros)
#  Covariance grid 4 x 3 = 12 (topology x overlap)
#  Proportions     H*=1 (balanced); H*=0.5; H*=0.1
#  Algorithms      DeconRNASeq (LSEI), CIBERSORT, L-BFGS-B,
#                  Newton–Raphson, Marquardt–Levenberg
#                  (barycentre start; no NNLS)
#  Replicates (n)  50
#  Total scenarios 4 x 3 x 3 = 36
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Rscript scripts/fig03_variance_driven.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/fig03/hybrid_config.rds (design + graph generator settings)
#  output/fig03/hybrid_{config,descriptors,theta,benchmark}.rds
#  output/fig03/density_visualisations/{mean_signature,scenario_metrics,
#    pairwise_hellinger,loglik_surface_p,loglik_rgl,network_topologies,
#    latent_projections}.pdf
#  output/fig03/performance_visualisations/{heatmap_rmse,heatmap_aitchison,raincloud,forest,
#    similarity,solver_dots,runtime,memory}.pdf
#  output/fig03/ggplot_rds/*.rds
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
DENSITY_DIR <- file.path(OUT_DIR, "density_visualisations")
PERF_DIR <- file.path(OUT_DIR, "performance_visualisations")
GGPLOT_RDS_DIR <- file.path(OUT_DIR, "ggplot_rds")
dir.create(DENSITY_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PERF_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(GGPLOT_RDS_DIR, recursive = TRUE, showWarnings = FALSE)

POSTPROCESS_ONLY <- identical(
  Sys.getenv("FIG03_POSTPROCESS_ONLY", "0"),
  "1"
)

.ui_h1("Figure 03 · Covariance-driven scenario")

N_REPL <- 100L
SEED <- 20260807L
set.seed(SEED)
N_MC_OVERLAP <- 4000L

N_GENES <- 20L
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
  scale_free = "#4C72B0",
  stochastic_block_model = "#C44E52"
)
ct_nms <- names(celltype_palette)

GRAPH_MODELS <- c("scale_free", "stochastic_block_model")
# One SPD completion per topology (fixed spectral cushion), then a
# global covariance scale s so MixSim BarOmega matches the target.
# Structural zeros of Omega_j are preserved because Omega(s) = Omega/s.
PRECISION_SHIFT_BASE <- 0.2
PRECISION_SCALE <- 0.3
PROP_INHIBITORY <- 0.5
OVERLAP_TARGET <- c(low = 0.02, moderate = 0.08, high = 0.20)

# Explicit generator settings (also written onto hybrid_config.rds).
GRAPH_PARAM_LIBRARY <- list(
  scale_free = list(
    power = 1,
    edges_per_node = 1L
  ),
  stochastic_block_model = list(
    block_prob = c(0.5, 0.25, 0.25),
    p_within = 0.25,
    p_between = 0.01
  )
)
GRAPH_CT3 <- "scale_free"

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
#   Fixed mean signature (exact Gram). Factorial graphs x target BarOmega.
# ==============================================================================

if (isTRUE(POSTPROCESS_ONLY)) {
  .ui_info("Post-process only: skipping generative model and ADEMP refit.")
} else {
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
    graph_ct2 = GRAPH_MODELS
  )
  topology_grid$graph_ct3 <- GRAPH_CT3
  overlap_grid <- tibble::tibble(
    overlap_label = names(OVERLAP_TARGET),
    overlap_target = unname(OVERLAP_TARGET)
  )

  cov_grid <- tidyr::expand_grid(topology_grid, overlap_grid)
  .ui_info(
    "Covariance grid: {.val {nrow(topology_grid)}} topologies x {.val {nrow(overlap_grid)}} overlap targets = {.val {nrow(cov_grid)}}."
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

  # ==============================================================================
  # SECTION 1a · Adjacency skeletons ----
  #   One undirected graph per cell type on each 2-by-2 topology row.
  # ==============================================================================

  .ui_info(
    "Drawing one adjacency triple per topology ({.val {nrow(topology_grid)}} graphs)."
  )
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
            graph_params = GRAPH_PARAM_LIBRARY[[models[[j]]]]
          )
        )
      })
    }
  )

  # ==============================================================================
  # SECTION 1b · Signed precision completion ----
  #   simulate_hierarchical_grn_moments() fills Omega_j on the skeleton,
  #   then inverts to Sigma_j at a fixed spectral cushion.
  # ==============================================================================

  n_topo <- nrow(topology_grid)
  .ui_info(
    "Completing signed precisions for each topology ({.val {n_topo}})."
  )
  base_draws <- lapply(seq_len(n_topo), function(topo_idx) {
    .ui_info("Base precision {.val {topo_idx}}/{.val {n_topo}}.")
    row <- topology_grid[topo_idx, , drop = FALSE]
    moments <- withr::with_seed(
      SEED + 100L * topo_idx,
      simulate_hierarchical_grn_moments(
        n_genes = N_GENES,
        n_celltypes = N_CELLTYPES,
        mean_scale = MEAN_SCALE,
        target_gram = TARGET_GRAM,
        precision_shift = PRECISION_SHIFT_BASE,
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
      ),
      graph_ct1 = row$graph_ct1,
      graph_ct2 = row$graph_ct2,
      graph_ct3 = row$graph_ct3
    )
  })

  # ==============================================================================
  # SECTION 1c · Scale covariances to target BarOmega ----
  #   Sigma_j |-> s Sigma_j (zeros of the completed Omega_j kept).
  # ==============================================================================

  n_cov <- nrow(cov_grid)
  .ui_info(
    "Scaling each topology to low / moderate / high MixSim BarOmega."
  )
  cov_draws <- lapply(seq_len(n_cov), function(i) {
    if (i == 1L || i %% 4L == 0L || i == n_cov) {
      .ui_info("Overlap scale {.val {i}}/{.val {n_cov}}.")
    }
    row <- cov_grid[i, , drop = FALSE]
    topo_idx <- which(
      topology_grid$graph_ct1 == row$graph_ct1 &
        topology_grid$graph_ct2 == row$graph_ct2 &
        topology_grid$graph_ct3 == row$graph_ct3
    )[[1L]]
    base <- base_draws[[topo_idx]]
    calibrated <- scale_covariances_to_overlap(
      mu = mu,
      sigma = base$sigma,
      p = p_balanced,
      target = row$overlap_target,
      n_mc = N_MC_OVERLAP,
      seed = SEED + i,
      verbose = FALSE
    )
    list(
      sigma = calibrated$sigma,
      theta = calibrated$Theta,
      adjacency = base$adjacency,
      covariance_scale = calibrated$scale,
      baromega = calibrated$baromega,
      overlap_target = calibrated$target
    )
  })

  # ==============================================================================
  # SECTION 1d · Scenario table and descriptors ----
  #   Cross 12 covariance draws with three Shannon compositions (36 rows).
  # ==============================================================================

  .ui_info("Building scenario descriptors (MixSim, f_cov, AIRM).")
  scenario_config_3 <- purrr::map_dfr(
    seq_len(n_cov),
    function(i) {
      if (i == 1L || i %% 8L == 0L || i == n_cov) {
        .ui_info("Scenario descriptors {.val {i}}/{.val {n_cov}}.")
      }
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
          adjacency = draw$adjacency,
          include_mixsim = TRUE
        )
        tibble::tibble(
          n_genes = N_GENES,
          n_celltypes = N_CELLTYPES,
          proportion_name = prop_name,
          entropy = round(compute_shannon_entropy(p), 3),
          graph_ct1 = row$graph_ct1,
          graph_ct2 = row$graph_ct2,
          graph_ct3 = row$graph_ct3,
          graph_generators = list(list(
            celltype_1 = list(
              graph_model = row$graph_ct1,
              graph_params = GRAPH_PARAM_LIBRARY[[row$graph_ct1]]
            ),
            celltype_2 = list(
              graph_model = row$graph_ct2,
              graph_params = GRAPH_PARAM_LIBRARY[[row$graph_ct2]]
            ),
            celltype_3 = list(
              graph_model = row$graph_ct3,
              graph_params = GRAPH_PARAM_LIBRARY[[row$graph_ct3]]
            )
          )),
          overlap_label = row$overlap_label,
          overlap_target = row$overlap_target,
          covariance_scale = draw$covariance_scale,
          mixsim_baromega = described$descriptors$mixsim_baromega,
          f_cov = described$descriptors$f_cov,
          kappa_sigma_p = described$descriptors$kappa_sigma_p,
          riemannian_sigma = described$descriptors$riemannian_sigma,
          true_theta = list(list(
            p = p,
            mu = mu,
            sigma = draw$sigma,
            Theta = draw$theta,
            adjacency = draw$adjacency,
            graph_generators = list(
              celltype_1 = list(
                graph_model = row$graph_ct1,
                graph_params = GRAPH_PARAM_LIBRARY[[row$graph_ct1]]
              ),
              celltype_2 = list(
                graph_model = row$graph_ct2,
                graph_params = GRAPH_PARAM_LIBRARY[[row$graph_ct2]]
              ),
              celltype_3 = list(
                graph_model = row$graph_ct3,
                graph_params = GRAPH_PARAM_LIBRARY[[row$graph_ct3]]
              )
            )
          ))
        )
      })
    }
  )
  .ui_success(
    "Config built: {.val {nrow(scenario_config_3)}} scenarios."
  )

  # ==============================================================================
  # SECTION 1e · Write design table ----
  # ==============================================================================

  overlap_check <- scenario_config_3[
    c(
      "proportion_name",
      "overlap_label",
      "overlap_target",
      "mixsim_baromega",
      "f_cov",
      "covariance_scale",
      "riemannian_sigma"
    )
  ]
  .ui_info("Realised BarOmega / f_cov by topology and overlap target:")
  print(overlap_check)
  scenario_config_3$ID <- paste0(
    "V",
    seq_len(nrow(scenario_config_3))
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
  write_simulation_artefacts(
    hybrid_out,
    dir = OUT_DIR,
    stem = "hybrid",
    config = scenario_config_3
  )
}

# ==============================================================================
# SECTION 3 · VISUALISATIONS ----
#   Folders match fig02: density_visualisations, performance_visualisations,
#   ggplot_rds. Nine pages (composition x overlap) with a 2-by-2 of CT1/CT2
#   topologies on each performance page.
# ==============================================================================

artefacts <- read_simulation_artefacts(
  OUT_DIR,
  "hybrid",
  assemble = TRUE
)
artefacts$theta <- readRDS(file.path(OUT_DIR, "hybrid_theta.rds"))
artefacts$config <- .relevel_scenario_table(artefacts$config)

old_root_pdfs <- file.path(
  OUT_DIR,
  c(
    "fig03_raincloud.pdf",
    "fig03_forest.pdf",
    "fig03_metric_dots.pdf",
    "fig03_network_topologies.png"
  )
)
unlink(old_root_pdfs[file.exists(old_root_pdfs)])

theta_one <- .unwrap_true_theta(artefacts$theta$true_theta[[1L]])

.ui_info("Drawing mean-signature heatmap.")
save_mean_signature_heatmap(
  theta_one,
  file.path(DENSITY_DIR, "mean_signature.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
.ui_success("Saved {.file mean_signature.pdf}.")

if (requireNamespace("funkyheatmap", quietly = TRUE)) {
  .ui_info("Drawing scenario-metric funkyheatmap.")
  save_scenario_metrics_funkyheatmap(
    artefacts,
    file.path(DENSITY_DIR, "scenario_metrics.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file scenario_metrics.pdf}.")
} else {
  .ui_warn("{.pkg funkyheatmap} not available; skipping scenario metrics.")
}

.ui_info("Drawing pairwise Hellinger 3-by-3 heatmap.")
save_pairwise_network_distance_heatmap(
  artefacts,
  file.path(DENSITY_DIR, "pairwise_hellinger.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
.ui_success("Saved {.file pairwise_hellinger.pdf}.")

if (requireNamespace("gridExtra", quietly = TRUE)) {
  .ui_info("Drawing ALR log-likelihood surfaces (9 pages, 2-by-2).")
  save_bivariate_loglik_surface_p_book(
    artefacts$config,
    artefacts$theta,
    file.path(DENSITY_DIR, "loglik_surface_p.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file loglik_surface_p.pdf}.")
  save_bivariate_expected_loglik_surface_p_book(
    artefacts$config,
    artefacts$theta,
    file.path(DENSITY_DIR, "loglik_expected_surface_p.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file loglik_expected_surface_p.pdf}.")
  if (
    requireNamespace("rgl", quietly = TRUE) &&
      requireNamespace("png", quietly = TRUE)
  ) {
    .ui_info("Drawing ALR log-likelihood rgl surfaces (9 pages).")
    save_bivariate_loglik_rgl_book(
      artefacts$config,
      artefacts$theta,
      file.path(DENSITY_DIR, "loglik_rgl.pdf"),
      html_file = file.path(DENSITY_DIR, "loglik_rgl.html")
    )
    .ui_success("Saved {.file loglik_rgl.pdf}.")
  } else {
    .ui_warn("{.pkg rgl} / {.pkg png} not available; skipping ALR rgl.")
  }
} else {
  .ui_warn("{.pkg gridExtra} not available; skipping ALR log-likelihood.")
}

if (requireNamespace("igraph", quietly = TRUE)) {
  .ui_info("Drawing 3-page network topology book.")
  vignette_network_png <- file.path(
    "vignettes",
    "figures",
    "fig_network_topologies.png"
  )
  save_network_topology_book(
    artefacts,
    file.path(DENSITY_DIR, "network_topologies.pdf"),
    data_rds = GGPLOT_RDS_DIR,
    png_file = file.path(DENSITY_DIR, "network_topologies.png")
  )
  file.copy(
    file.path(DENSITY_DIR, "network_topologies.png"),
    vignette_network_png,
    overwrite = TRUE
  )
  .ui_success("Saved {.file network_topologies.pdf}.")
} else {
  .ui_warn("{.pkg igraph} not available; skipping network book.")
}

.ui_info("Drawing RMSE and Aitchison tile heatmaps.")
save_bivariate_metric_heatmaps(
  artefacts,
  PERF_DIR,
  data_rds = GGPLOT_RDS_DIR
)
.ui_success("Saved RMSE heatmap PDF.")

if (requireNamespace("ggdist", quietly = TRUE)) {
  .ui_info("Drawing 9-page raincloud book.")
  save_bivariate_raincloud_book(
    artefacts,
    file.path(PERF_DIR, "raincloud.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file raincloud.pdf}.")
  .ui_info("Drawing 9-page quantile-dot raincloud book.")
  save_bivariate_dotsinterval_book(
    artefacts,
    file.path(PERF_DIR, "dotsinterval.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file dotsinterval.pdf}.")
} else {
  .ui_warn("{.pkg ggdist} not available; skipping raincloud book.")
}

.ui_info("Drawing 9-page Wald forest book.")
save_bivariate_forest_book(
  artefacts,
  file.path(PERF_DIR, "forest.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
.ui_success("Saved {.file forest.pdf}.")

if (requireNamespace("ggdist", quietly = TRUE)) {
  .ui_info("Drawing KKT raincloud book.")
  save_bivariate_kkt_raincloud_book(
    artefacts,
    file.path(PERF_DIR, "kkt_raincloud.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file kkt_raincloud.pdf}.")
}

if (requireNamespace("cowplot", quietly = TRUE)) {
  .ui_info("Drawing stacked convergence book.")
  save_bivariate_convergence_book(
    artefacts,
    file.path(PERF_DIR, "convergence_stacked.pdf"),
    data_rds = GGPLOT_RDS_DIR,
    icon_dir = "temp_logos"
  )
  .ui_success("Saved {.file convergence_stacked.pdf}.")
}

.ui_info("Drawing Q-Q and KS normality books.")
save_bivariate_qq_book(
  artefacts,
  file.path(PERF_DIR, "qq_normal.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
save_bivariate_ks_book(
  artefacts,
  file.path(PERF_DIR, "ks_normal.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
save_bivariate_chi2_book(
  artefacts,
  file.path(PERF_DIR, "qq_chi2.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
.ui_success(
  "Saved {.file qq_normal.pdf}, {.file ks_normal.pdf} and {.file qq_chi2.pdf}."
)

if (length(unique(artefacts$monte_carlo$algorithm)) >= 2L) {
  if (
    requireNamespace("ggdendro", quietly = TRUE) &&
      requireNamespace("cowplot", quietly = TRUE) &&
      requireNamespace("gridExtra", quietly = TRUE)
  ) {
    .ui_info("Drawing 9-page similarity book with dendrograms.")
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

.ui_info("Drawing 9-page RMSE / Aitchison solver-dot book.")
save_bivariate_solver_dots_book(
  artefacts,
  file.path(PERF_DIR, "solver_dots.pdf"),
  data_rds = GGPLOT_RDS_DIR
)
.ui_success("Saved {.file solver_dots.pdf}.")

if (
  requireNamespace("EMMIXmfa", quietly = TRUE) &&
    requireNamespace("mclust", quietly = TRUE) &&
    requireNamespace("cowplot", quietly = TRUE)
) {
  .ui_info("Drawing 16-page 2D latent-projection book.")
  save_hybrid_latent_projection_book(
    artefacts,
    file.path(DENSITY_DIR, "latent_projections.pdf"),
    n = 300L,
    seed = SEED,
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file latent_projections.pdf}.")
} else {
  .ui_warn(
    "Latent-projection book skipped (need EMMIXmfa, mclust, cowplot)."
  )
}

if (requireNamespace("ggdist", quietly = TRUE)) {
  .ui_info("Drawing solver runtime rainclouds (3 pages).")
  save_hybrid_runtime_book(
    artefacts,
    file.path(PERF_DIR, "runtime.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file runtime.pdf}.")
  .ui_info("Drawing solver memory rainclouds (3 pages).")
  save_hybrid_memory_book(
    artefacts,
    file.path(PERF_DIR, "memory.pdf"),
    data_rds = GGPLOT_RDS_DIR
  )
  .ui_success("Saved {.file memory.pdf}.")
} else {
  .ui_warn("{.pkg ggdist} not available; skipping runtime / memory books.")
}

.ui_success(
  "Done. Outputs in {.path {normalizePath(OUT_DIR, mustWork = FALSE)}}."
)
