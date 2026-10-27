###############################################################################
###                                                                         ###
###     FIGURE 03 - COVARIANCE-DRIVEN SCENARIO (article Section 2.3)        ###
###     J = 3 cell types, G = 20 genes, network-constrained covariances     ###
###                                                                         ###
###############################################################################
#
# Run from the repository root; stdout and stderr go to logs/:
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore \
#     scripts/fig03_variance_covariance_driven.R \
#     > "logs/fig03_$(date +%F)_covariance_driven.log" 2>&1 &
#
# Run size (environment variables; see
# scripts/benchmark_metrics/benchmark_artefacts.R):
#   DECOVART_N_REPLICATES  Monte Carlo bulks per scenario   (default 100)
#   DECOVART_N_BOOT        bootstrap refits per bulk          (default 200)
#   DECOVART_N_REFERENCE   individuals of the reference panel (default N)
#   DECOVART_CORES         workers for the per-sample loop    (default 1)
#                          (workers load the installed DeCovarT, so run
#                          devtools::install() before a parallel run)
#   DECOVART_POSTPROCESS_ONLY=1 redraws the figures from the saved benchmark.
#
# Design (36 scenarios, IDs V1 to V36)
#  Mean signature  one Gram: cosines 0.9 (types 1-2) and 0.1 (type 3),
#                  mean scale 10, non-negative
#  Topology        type 1 x type 2 in {scale-free, stochastic block model};
#                  type 3 scale-free
#  Overlap         average MixSim overlap 0.02, 0.08 or 0.20 at the balanced
#                  composition, by a global covariance scale
#  Composition     normalised Shannon entropy 1, 0.5 or 0.1
#
# Methods
#  lsei, NNLS, RLM,       first-generation baselines on the true means: lsei
#  CIBERSORT              keeps the simplex at every iteration, NNLS fits a
#                         free scale then renormalises, RLM is the robust
#                         linear model of Monaco et al.
#  Newton-Raphson,        DeCovarT on the true moments; intervals: Wald, ILR,
#  L-BFGS-B, Firth        profile, bootstrap (+ oracle for Newton-Raphson)
#  MEAD, RNA-Sieve,       authors' code (scripts/benchmark_metrics/
#  ReDeconv               second_generation_methods.R); MEAD and RNA-Sieve
#                         estimate the reference from an independent panel of
#                         R purified profiles per type, fitted on the cohort
#                         of the N bulks of a scenario; ReDeconv takes the true
#                         means and standard deviations
#  DECALS                 no public code: not run (its point estimate is lsei)
#
# Outputs (output/fig03/)
#  fig03_benchmark.rds, fig03_plot_data.rds,
#  density_visualisations/*.pdf, performance_visualisations/*.pdf
###############################################################################

devtools::load_all(".", quiet = TRUE)
source(file.path("scripts", "benchmark_metrics", "benchmark_artefacts.R"))
source(file.path("scripts", "benchmark_metrics", "second_generation_methods.R"))
source(file.path("scripts", "visualisation", "load_visualisation.R"))
stopifnot(
  requireNamespace("igraph", quietly = TRUE),
  requireNamespace("e1071", quietly = TRUE)
)

OUT_DIR <- file.path("output", "fig03")
BENCHMARK_FILE <- file.path(OUT_DIR, "fig03_benchmark.rds")
SEED <- 20260807L
CORES <- as.integer(Sys.getenv("DECOVART_CORES", "1"))
settings <- benchmark_settings(default_replicates = 100L, default_boot = 200L)

N_MC_OVERLAP <- 4000L
N_GENES <- 20L
N_CELLTYPES <- 3L
MEAN_SCALE <- 10
ITMAX <- 200L
EPSILON <- 1e-4
ct_nms <- paste0("celltype_", seq_len(N_CELLTYPES))

GRAPH_MODELS <- c("scale_free", "stochastic_block_model")
# One positive-definite completion per topology (fixed spectral cushion),
# then a global covariance scale so the MixSim overlap hits its target;
# the zeros of every precision matrix are kept.
PRECISION_SHIFT_BASE <- 0.2
PRECISION_SCALE <- 0.3
PROP_INHIBITORY <- 0.5
OVERLAP_TARGET <- c(low = 0.02, moderate = 0.08, high = 0.20)
GRAPH_PARAM_LIBRARY <- list(
  scale_free = list(power = 1, edges_per_node = 1L),
  # Three contiguous blocks of near-equal size (7, 6, 7 genes), the same
  # rule as fig05 and as the network plot.
  stochastic_block_model = list(
    block_sizes = tabulate(cut(seq_len(N_GENES), 3L, labels = FALSE)),
    p_within = 0.25,
    p_between = 0.01
  )
)
GRAPH_CT3 <- "scale_free"
TARGET_GRAM <- matrix(
  c(1, 0.9, 0.1, 0.9, 1, 0.1, 0.1, 0.1, 1),
  nrow = 3L,
  dimnames = list(ct_nms, ct_nms)
)

.ui_h1("Figure 03: covariance-driven scenario")

if (!settings$postprocess_only) {
  set.seed(SEED)
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
          proportions = prop_name,
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

  scenario_config_3$ID <- paste0("V", seq_len(nrow(scenario_config_3)))

  interval_menu <- function(methods, optimiser = "newton") {
    list(methods = methods, n_boot = settings$n_boot, optimiser = optimiser)
  }
  solver <- list(epsilon = EPSILON, itmax = ITMAX, initial_p = "barycentre")
  deconvolution_functions_3 <- list(
    "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
    "NNLS" = list(FUN = deconvolute_ratios_nnls),
    "RLM" = list(FUN = deconvolute_ratios_rlm),
    "CIBERSORT" = list(FUN = deconvolute_ratios_cibersort),
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

  .ui_info(
    "{.val {nrow(scenario_config_3)}} scenarios, N = {.val {settings$n_replicates}}, B = {.val {settings$n_boot}}, R = {.val {settings$n_reference}}."
  )
  started <- proc.time()[["elapsed"]]
  hybrid_out <- run_simulation_benchmark(
    scenario_config = scenario_config_3,
    deconvolution_functions = deconvolution_functions_3,
    n = settings$n_replicates,
    cores = CORES,
    verbose = TRUE,
    return_data = TRUE,
    reference = list(
      n_individuals = settings$n_reference,
      n_cells = 1L,
      seed = SEED + 5000L
    )
  )
  .ui_info("Second-generation methods on the same bulks.")
  competitors <- run_second_generation_batch(
    competitor_scenarios(hybrid_out, scenario_config_3)
  )
  hybrid_out <- score_competitors(
    hybrid_out,
    competitors,
    scenario_config_3,
    model = "population"
  )
  .ui_success(
    "Benchmark finished in {.val {round(proc.time()[['elapsed']] - started)}} s."
  )
  save_benchmark(
    hybrid_out,
    scenario_config_3,
    benchmark_meta(SEED, settings),
    BENCHMARK_FILE
  )
}

benchmark <- read_benchmark(BENCHMARK_FILE)
plot_data <- figure_plot_data(benchmark, figure = "fig03")
saveRDS(plot_data, file.path(OUT_DIR, "fig03_plot_data.rds"))
render_figure_books(plot_data, OUT_DIR)
.ui_success("Figure 03 written to {.file {OUT_DIR}}.")
