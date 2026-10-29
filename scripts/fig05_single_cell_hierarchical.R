###############################################################################
###                                                                         ###
###     FIGURE 05 - SINGLE-CELL REFERENCE, TWO LINEAGES (article 2.4)       ###
###     J = 4 cell types in two lineages, bulk = sum of n_i cells           ###
###                                                                         ###
###############################################################################
#
# Run from the repository root; stdout and stderr go to logs/. OpenBLAS
# 0.3.20 needs OPENBLAS_CORETYPE on this machine (see the check below):
#
#   mkdir -p logs
#   export OPENBLAS_CORETYPE=SkylakeX OPENBLAS_NUM_THREADS=1
#   nohup Rscript --no-save --no-restore \
#     scripts/fig05_single_cell_hierarchical.R \
#     > "logs/fig05_$(date +%F)_single_cell.log" 2>&1 &
#
# Run size, hard-coded below: N_REPLICATES Monte Carlo bulks per scenario,
# N_REFERENCE individuals of the reference panel (MEAD and RNA-Sieve only),
# CORES workers for the N bulks of each scenario (they load the installed
# DeCovarT, so run devtools::install() first). POSTPROCESS_ONLY <- TRUE redraws the
# figures from the saved benchmark.
#
# Generative model: the single-cell mixture of the article,
#   y = a sum_j T_j,  T_j ~ N(m_j S_j theta_j, m_j S_j^2 Sigma^theta_j),
#   m ~ Multinomial(n_i, p); the target is the realised fraction m / n_i.
#
# Design (64 crossed scenarios: 2^5 at a = 5 and the same 2^5 at a = 1)
#  Lineages        types 1-2 and types 3-4; mean cosines 0.8 within a
#                  lineage and 0.3 between, non-negative means with no
#                  gene expressed in one type only
#  Networks        one skeleton per lineage, shared exactly by its two
#                  siblings: a stochastic block model with three equal
#                  blocks for types 1-2, a scale-free graph for types 3-4;
#                  siblings share the lineage's signed weights except for
#                  5% of edges whose sign each flips independently
#  Composition     balanced (1/4 each) or unbalanced within each lineage
#                  (0.4, 0.1, 0.4, 0.1)
#  Overlap         mean MixSim overlap of the two sibling pairs 0.01 or 0.08
#  Cells per bulk  n_i = 1e3 or 1e5
#  Cell sizes      equal, or lineage 3-4 five times larger (S = 1, 1, 5, 5);
#                  the size scales the mean and, squared, the covariance
#  Genes           G = 20 or 200 (the bias of phi falls as (G - J) / G)
#  Platform scale  a = 1 or 5; a = 1 isolates the scale assumption of
#                  RNA-Sieve (p-hat of DeCovarT does not depend on a), and
#                  the crossed a = 1 half checks that invariance
#                  numerically
#
# Methods
#  lsei                   simplex least squares with the bulk scale fixed at
#                         1, misspecified here since the scale is c = a n
#  NNLS                   least squares with a free bulk scale, renormalised
#                         (the start of DeCovarT and RNA-Sieve)
#  RLM                    robust linear model (MASS::rlm), renormalised
#  CIBERSORT              nu-SVR on the true means
#  Newton-Raphson,        DeCovarT single-cell fit (cell fractions, known
#  L-BFGS-B, Firth        cell sizes, raw covariances); intervals: Wald, ILR,
#                         profile (ML and REML), oracle
#  MEAD, RNA-Sieve,       authors' code; MEAD and RNA-Sieve read a reference
#  ReDeconv               panel of R individuals x 20 cells per type drawn
#                         independently of the bulks; ReDeconv takes the true
#                         means and standard deviations
#
# Outputs (output/fig05/)
#  fig05_benchmark.rds, fig05_plot_data.rds,
#  density_visualisations/*.pdf, performance_visualisations/*.pdf
###############################################################################

# Run size and output.
N_REPLICATES <- 200L
N_REFERENCE <- 200L
# The three figure scripts run at the same time on half the machine:
# a third of half the cores each (14 on 88 logical cores).
CORES <- (parallel::detectCores() %/% 2L) %/% 3L
POSTPROCESS_ONLY <- FALSE

# Locally installed build, not devtools::load_all() of this tree.
library(DeCovarT)
# A faulty BLAS kernel returns wrong Cholesky factors and inverses
# (OpenBLAS 0.3.20 on AVX-512 / AMX CPUs).
local({
  a <- 0.5^abs(outer(seq_len(64L), seq_len(64L), "-"))
  if (
    max(abs(crossprod(chol(a)) - a)) > 1e-10 ||
      max(abs(solve(a) %*% a - diag(64L))) > 1e-8
  ) {
    stop(
      "BLAS returns wrong chol()/solve(); set OPENBLAS_CORETYPE=SkylakeX ",
      "before starting R.",
      call. = FALSE
    )
  }
})
source(file.path("scripts", "benchmark_metrics", "benchmark_artefacts.R"))
source(file.path("scripts", "benchmark_metrics", "second_generation_methods.R"))
source(file.path("scripts", "visualisation", "load_visualisation.R"))
stopifnot(
  requireNamespace("igraph", quietly = TRUE),
  requireNamespace("e1071", quietly = TRUE)
)

OUT_DIR <- file.path("output", "fig05")
BENCHMARK_FILE <- file.path(OUT_DIR, "fig05_benchmark.rds")
SEED <- 20261024L

N_CELLTYPES <- 4L
ct_nms <- paste0("celltype_", seq_len(N_CELLTYPES))
SIBLINGS <- rbind(c(1L, 2L), c(3L, 4L))
MEAN_SCALE <- 10
PRECISION_SHIFT <- 0.2
PRECISION_SCALE <- 0.3
PROP_INHIBITORY <- 0.5
N_MC_OVERLAP <- 4000L
N_CELLS_REFERENCE <- 20L
# One skeleton per lineage: a stochastic block model with three blocks of
# equal size for types 1-2, a scale-free graph for types 3-4.
LINEAGE <- c(1L, 1L, 2L, 2L)
LINEAGE_GRAPH <- c("stochastic_block_model", "scale_free")
# Siblings keep the signed weights of their lineage except on a random 5%
# of edges, whose sign each sibling flips independently. With independent
# signs (or 20% flips) the sibling covariances differ so much at G = 200 that
# a sibling overlap of 0.08 is out of reach (maximum 0.013, or 0.030); with
# 5% every target is met inside the search range of the covariance scale.
SIGN_FLIP <- 0.05
graph_params <- function(model, n_genes) {
  switch(
    model,
    stochastic_block_model = list(
      block_sizes = tabulate(cut(seq_len(n_genes), 3L, labels = FALSE)),
      p_within = 0.25,
      p_between = 0.01
    ),
    scale_free = list(power = 1, edges_per_node = 1L)
  )
}
TARGET_GRAM <- matrix(0.3, N_CELLTYPES, N_CELLTYPES)
TARGET_GRAM[1:2, 1:2] <- 0.8
TARGET_GRAM[3:4, 3:4] <- 0.8
diag(TARGET_GRAM) <- 1
dimnames(TARGET_GRAM) <- list(ct_nms, ct_nms)

COMPOSITIONS <- list(
  balanced = rep(0.25, N_CELLTYPES),
  unbalanced = c(0.4, 0.1, 0.4, 0.1)
)
OVERLAP_TARGET <- c(low = 0.01, high = 0.08)
N_CELLS <- c(1e3, 1e5)
CELL_SIZES <- list(equal = c(1, 1, 1, 1), lineage = c(1, 1, 5, 5))
N_GENES <- c(20L, 200L)
# Both platform scales on the full factorial grid (2^5 x 2 = 64 scenarios): the
# a = 1 half is redundant with a = 5 up to the nuisance scale, kept as a
# numerical check.
PLATFORM_SCALE <- c(1, 5)
ITMAX <- 100L
EPSILON <- 1e-8

# Precisions of the four types: the two siblings of a lineage share its
# skeleton exactly and differ only by their edge weights.
lineage_moments <- function(n_genes, seed) {
  withr::with_seed(seed, {
    skeletons <- lapply(LINEAGE_GRAPH, function(model) {
      generate_random_network_skeleton(
        n_genes = n_genes,
        graph_model = model,
        graph_params = graph_params(model, n_genes)
      )
    })
    adjacency <- skeletons[LINEAGE]
    lineage_weights <- lapply(skeletons, function(skeleton) {
      assign_iid_signed_weights(skeleton, PROP_INHIBITORY, PRECISION_SCALE)
    })
    precision <- array(
      0,
      dim = c(n_genes, n_genes, N_CELLTYPES),
      dimnames = list(paste0("gene_", seq_len(n_genes)), NULL, ct_nms)
    )
    dimnames(precision)[[2L]] <- dimnames(precision)[[1L]]
    for (j in seq_len(N_CELLTYPES)) {
      weights <- lineage_weights[[LINEAGE[[j]]]]
      edges <- which(upper.tri(weights) & weights != 0)
      flip <- edges[stats::runif(length(edges)) < SIGN_FLIP]
      weights[flip] <- -weights[flip]
      weights[lower.tri(weights)] <- t(weights)[lower.tri(weights)]
      precision[,, j] <- build_normalised_precision(weights, PRECISION_SHIFT)
    }
    mu <- generate_mean_signature_matrix(
      n_genes = n_genes,
      n_celltypes = N_CELLTYPES,
      mean_scale = MEAN_SCALE,
      target_gram = TARGET_GRAM,
      seed = seed,
      nonnegative = TRUE,
      celltype_names = ct_nms
    )
    list(
      mean_profiles = mu,
      covariance_matrices = build_covariance_array_from_precision(precision),
      adjacency = adjacency
    )
  })
}

# One reference per gene count, overlap target and cell-size vector. The
# covariance scale is calibrated once per (G, overlap) on equal sizes:
# a common size multiplies the means by S and the covariances by S^2, so the
# sibling overlap does not depend on the sizes.
build_single_cell_config <- function() {
  references <- list()
  for (g in seq_along(N_GENES)) {
    moments <- lineage_moments(N_GENES[[g]], SEED + 100L * g)
    unit <- rescale_single_cell_moments(moments, cell_size = 1)
    for (o in seq_along(OVERLAP_TARGET)) {
      calibrated <- scale_covariances_to_overlap(
        mu = unit$generative_mean,
        sigma = unit$generative_sigma,
        p = COMPOSITIONS$balanced,
        target = OVERLAP_TARGET[[o]],
        n_mc = N_MC_OVERLAP,
        seed = SEED + 10L * g + o,
        pairs = SIBLINGS
      )
      scaled <- moments
      scaled$covariance_matrices <- moments$covariance_matrices *
        calibrated$scale
      for (size in names(CELL_SIZES)) {
        references[[length(references) + 1L]] <- list(
          n_genes = N_GENES[[g]],
          overlap_label = names(OVERLAP_TARGET)[[o]],
          overlap_target = OVERLAP_TARGET[[o]],
          sibling_overlap = calibrated$baromega,
          covariance_scale = calibrated$scale,
          cell_size_label = size,
          reference = rescale_single_cell_moments(scaled, CELL_SIZES[[size]]),
          adjacency = moments$adjacency
        )
      }
    }
  }
  design <- tidyr::expand_grid(
    reference_idx = seq_along(references),
    proportions = names(COMPOSITIONS),
    n_cells = N_CELLS,
    platform_scale = PLATFORM_SCALE
  )
  purrr::pmap_dfr(
    design,
    function(reference_idx, proportions, n_cells, platform_scale) {
      r <- references[[reference_idx]]
      ref <- r$reference
      p <- stats::setNames(COMPOSITIONS[[proportions]], ct_nms)
      tibble::tibble(
        n_genes = r$n_genes,
        n_celltypes = N_CELLTYPES,
        proportions = proportions,
        entropy = round(compute_shannon_entropy(p), 3),
        overlap_label = r$overlap_label,
        overlap_target = r$overlap_target,
        sibling_overlap = r$sibling_overlap,
        covariance_scale = r$covariance_scale,
        cell_size_label = r$cell_size_label,
        cell_size_ratio = max(ref$cell_size) / min(ref$cell_size),
        n_cells = n_cells,
        platform_scale = platform_scale,
        true_theta = list(list(
          p = p,
          mu = ref$mean_signature_matrix,
          sigma = ref$Sigma,
          generative_mean = ref$generative_mean,
          generative_sigma = ref$generative_sigma,
          n_cells = n_cells,
          platform_scale = platform_scale,
          cell_size = ref$cell_size,
          adjacency = r$adjacency
        ))
      )
    }
  ) |>
    dplyr::mutate(
      ID = sprintf(
        "S%d_%s_%s_n%s_%s_G%d_a%g",
        dplyr::row_number(),
        ifelse(.data$proportions == "balanced", "Ba", "Un"),
        ifelse(.data$overlap_label == "low", "O1", "O8"),
        ifelse(.data$n_cells < 1e4, "1e3", "1e5"),
        ifelse(.data$cell_size_label == "equal", "Se", "Sl"),
        .data$n_genes,
        .data$platform_scale
      )
    )
}

cli::cli_h1("Figure 05: single-cell reference, two lineages")

if (!POSTPROCESS_ONLY) {
  set.seed(SEED)
  scenario_config_5 <- build_single_cell_config()
  interval_menu <- function(methods, optimiser = "newton") {
    list(methods = methods, optimiser = optimiser)
  }
  menu <- c(
    "wald",
    "ilr",
    "ilr_reml",
    "profile",
    "profile_reml"
  )
  solver <- list(epsilon = EPSILON, itmax = ITMAX)
  deconvolution_functions_5 <- list(
    "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
    "NNLS" = list(FUN = deconvolute_ratios_nnls),
    "RLM" = list(FUN = deconvolute_ratios_rlm),
    "CIBERSORT" = list(FUN = deconvolute_ratios_cibersort),
    "Newton-Raphson" = list(
      FUN = deconvolute_single_cell,
      additional_parameters = c(solver, optimiser = "newton"),
      ci = interval_menu(c(menu, "oracle"))
    ),
    "L-BFGS-B" = list(
      FUN = deconvolute_single_cell,
      additional_parameters = c(solver, optimiser = "lbfgsb"),
      # Profile intervals refit from the bulk and do not depend on the
      # solver: they are reported once, under Newton-Raphson.
      ci = interval_menu(c("wald", "ilr", "ilr_reml"), "lbfgsb")
    ),
    "Firth" = list(
      FUN = fit_firth_single_cell,
      ci = interval_menu(c("wald", "ilr", "wald_reml", "ilr_reml"))
    )
  )
  cli::cli_alert_info(
    "{.val {nrow(scenario_config_5)}} scenarios, N = {.val {N_REPLICATES}}, R = {.val {N_REFERENCE}}, {.val {CORES}} workers."
  )
  # One pool of workers for the whole run: deconvolute_ratios() reuses a
  # plan that is already set, instead of starting CORES new R sessions
  # for every scenario (7 s each time for a pool of 44 workers).
  if (CORES > 1L) {
    future::plan(future::multisession, workers = CORES)
  }
  started <- proc.time()[["elapsed"]]
  single_cell_out <- run_simulation_benchmark(
    scenario_config = scenario_config_5,
    deconvolution_functions = deconvolution_functions_5,
    n = N_REPLICATES,
    cores = CORES,
    verbose = TRUE,
    likelihood = single_cell_likelihood(),
    return_data = TRUE,
    # Scenario i runs under seed SEED + i, so a failed scenario is rerun
    # alone with `scenarios = i`.
    seed = SEED,
    reference = list(
      n_individuals = N_REFERENCE,
      n_cells = N_CELLS_REFERENCE,
      seed = SEED + 5000L
    )
  )
  cli::cli_alert_info("Second-generation methods on the same bulks.")
  competitors <- run_second_generation_batch(
    competitor_scenarios(single_cell_out, scenario_config_5),
    workers = CORES
  )
  single_cell_out <- score_competitors(
    single_cell_out,
    competitors,
    scenario_config_5,
    likelihood = single_cell_likelihood()
  )
  cli::cli_alert_success(
    "Benchmark finished in {.val {round(proc.time()[['elapsed']] - started)}} s."
  )
  save_benchmark(
    single_cell_out,
    scenario_config_5,
    benchmark_meta(SEED, N_REPLICATES, N_REFERENCE),
    BENCHMARK_FILE
  )
}

benchmark <- readRDS(BENCHMARK_FILE)
plot_data <- figure_plot_data(benchmark, figure = "fig05")
saveRDS(plot_data, file.path(OUT_DIR, "fig05_plot_data.rds"))
render_figure_books(plot_data, OUT_DIR)
cli::cli_alert_success("Figure 05 written to {.file {OUT_DIR}}.")
