###############################################################################
###                                                                         ###
###     FIGURE 05 - BIVARIATE TOY, SINGLE-CELL MODEL (article Section 2.4)  ###
###       J = 2 cell types, G = 2 genes, known technical scales             ###
###                                                                         ###
###############################################################################
#
# Run from the repository root; stdout and stderr go to logs/ (see
# scripts/run_figure_benchmarks.sh for the OpenBLAS settings):
#
#   bash scripts/run_figure_benchmarks.sh fig05
#
# The single-cell twin of scripts/fig02_bivariate_toy.R. A bulk sums n_i
# cells, y_i ~ N(c_i mu p_i, phi_i S_W(p_i)), S_W(p) = sum_j p_j Sigma_W,j,
# with the platform scale fixed at a = 1, so c_i = phi_i = n_i. With G = 2
# genes a bulk cannot identify c_i or phi_i (the profiled likelihood has no
# maximum), so every scale-dependent method is given them, in two arms:
#  oracle   the true n_i of the bulk
#  nominal  the design value n; the bulks draw n_i log-normal around it
#           (coefficient of variation 30%), as a cell count that is guessed
#           rather than known.
#
# Factorial design (729 scenarios)
#  Proportions             balanced (1/2, 1/2); moderately unbalanced
#                          (17/20, 3/20); highly unbalanced (99/100, 1/100)
#  Cells per bulk n        100; 1,000; 10,000 (nominal)
#  Correlation, type 1     -0.8 to 0.8 by 0.2 (9 levels)
#  Correlation, type 2     -0.8 to 0.8 by 0.2 (9 levels)
#  Fixed                   per-cell means (20, 22) / (22, 20) (the small
#                          centroid gap of fig02); per-cell covariances
#                          Sigma_W,j = 500 R_j (unit-diagonal R_j). At
#                          n = 1,000 and the balanced composition, the bulk
#                          covariance per cell, S_W(p) / n, equals fig02's
#                          homoscedastic Sigma(p) = sum_j p_j^2 R_j.
#  Scenario ID             SB{index}_{Ba|Mo|Hi}_n{1e2|1e3|1e4}
#  Target                  the realised cell fraction m_i / n_i
#
# Methods (each of the first seven in the oracle and in the nominal arm)
#  lsei                  DeconRNASeq on y / c (simplex-constrained LS)
#  DWLS                  DWLS::solveDampenedWLS() on y / c; switched off
#                        (RUN_DWLS) because it fails at G = 2
#  Newton-Raphson        DeCovarT, trust-region Newton in ILR coordinates
#  L-BFGS-B              DeCovarT, L-BFGS-B on the closure of the orthant
#  Marquardt-Levenberg   DeCovarT, marqLevAlg in ILR coordinates
#  Firth                 DeCovarT with the Jeffreys penalty
#  DeCovarT diagonal     Newton-Raphson with the off-diagonals of Sigma_W,j
#                        set to zero: the model of RNA-Sieve at the true
#                        moments
#  NNLS                  non-negative least squares, free scale
#  RNA-Sieve, ReDeconv   the authors' Python code (corrected RNA-Sieve
#                        update, see scripts/benchmark_metrics/python/)
# RLM, CIBERSORT, MEAD and DECALS are not run: two genes give two
# observations, below what their fits or asymptotics need.
#
# Inference
#  Wald, ILR, oracle      every DeCovarT solver; profile under Newton-Raphson
#  Parametric bootstrap   Newton-Raphson, B = 199, on 8 scenarios:
#                         {Ba, Mo} x {n = 100, 1,000} x {(0, 0), (-0.8, 0.8)};
#                         the nominal arm also with n redrawn from its
#                         log-normal law in every replicate
#                         (`bootstrap_nprior`)
#  Boundary LRT           p_2 = 0 (chi-bar-square) under Newton-Raphson,
#                         both arms; size when the bulk holds no type-2
#                         cell, power otherwise
#
# Outputs (output/fig05/)
#  fig05_benchmark.rds   settings, design, truth, every metric table, the
#                        boundary LRTs
#  fig05_plot_data.rds   the data behind every figure page, keyed by stem
#  density_visualisations/, performance_visualisations/: one SVG and one
#  RDS (its plot data) per page
#
# Sourcing this file from a test or a vignette defines the builders only.
###############################################################################

# Per-cell variance: at n = 1,000 and p = (1/2, 1/2), S_W(p) / n equals
# fig02's sum_j p_j^2 R_j when 500 / 1000 * (R_1 + R_2) / 2 = (R_1 + R_2) / 4.
FIG05_CELL_VARIANCE <- 500
FIG05_N_CELLS_CV <- 0.3

#' Build the bivariate single-cell scenario configuration
#'
#' @return Tibble with one row per scenario, its descriptors and the
#'   `true_theta` list column read by [run_simulation_benchmark()].
build_bivariate_single_cell_config <- function(
  proportions = list(
    "balanced" = c(0.5, 0.5),
    "moderately unbalanced" = c(0.85, 0.15),
    "highly unbalanced" = c(0.99, 0.01)
  ),
  n_cells = c(100, 1000, 10000),
  corr_sequence = seq(-0.8, 0.8, 0.2),
  cell_variance = FIG05_CELL_VARIANCE,
  n_cells_cv = FIG05_N_CELLS_CV
) {
  proportion_list <- proportions
  genes <- c("gene_1", "gene_2")
  types <- c("celltype_1", "celltype_2")
  mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, types))
  design <- tidyr::expand_grid(
    proportions = names(proportions),
    n_cells = n_cells,
    correlation_celltype1 = corr_sequence,
    correlation_celltype2 = corr_sequence
  )
  p_code <- c(
    "balanced" = "Ba",
    "moderately unbalanced" = "Mo",
    "highly unbalanced" = "Hi"
  )
  purrr::pmap(
    design,
    function(
      proportions,
      n_cells,
      correlation_celltype1,
      correlation_celltype2
    ) {
      p <- stats::setNames(proportion_list[[proportions]], types)
      sigma <- array(0, c(2L, 2L, 2L), dimnames = list(genes, genes, types))
      for (j in 1:2) {
        rho <- c(correlation_celltype1, correlation_celltype2)[[j]]
        sigma[,, j] <- cell_variance * matrix(c(1, rho, rho, 1), 2L)
      }
      fisher <- fig05_fisher_split(p, mu, sigma, n_cells)
      # Overlap of the two pure bulks, N(mu_j, Sigma_W,j / n), per cell.
      overlap <- tryCatch(
        signif(
          MixSim::overlap(Pi = p, Mu = t(mu), S = sigma / n_cells)$BarOmega,
          3
        ),
        error = function(e) NA_real_
      )
      tibble::tibble(
        correlation_celltype1 = correlation_celltype1,
        correlation_celltype2 = correlation_celltype2,
        overlap = overlap,
        entropy = round(compute_shannon_entropy(p), 3),
        proportions = proportions,
        n_cells = n_cells,
        fisher_mean = fisher[["mean"]],
        fisher_cov = fisher[["cov"]],
        f_cov = fisher[["cov"]] / (fisher[["mean"]] + fisher[["cov"]]),
        true_theta = list(list(
          p = p,
          mu = mu,
          sigma = sigma,
          n_cells = n_cells,
          n_cells_cv = n_cells_cv,
          platform_scale = 1
        ))
      )
    }
  ) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      ID = sprintf(
        "SB%d_%s_n%s",
        dplyr::row_number(),
        p_code[.data$proportions],
        c("100" = "1e2", "1000" = "1e3", "10000" = "1e4")[
          as.character(.data$n_cells)
        ]
      ),
      .before = 1L
    )
}

#' Fisher information of p_1 at known scales, split into its mean and
#' covariance parts
#'
#' For y / c ~ N(mu p, S_W(p) / n) and p = (p_1, 1 - p_1), with
#' d = mu_1 - mu_2 and D = Sigma_W,1 - Sigma_W,2:
#' I_mean = n d' S_W(p)^-1 d and I_cov = tr[(S_W(p)^-1 D)^2] / 2. Only the
#' mean part grows with n.
fig05_fisher_split <- function(p, mu, sigma, n_cells) {
  s_w <- p[[1L]] * sigma[,, 1L] + p[[2L]] * sigma[,, 2L]
  d <- mu[, 1L] - mu[, 2L]
  a <- solve(s_w, sigma[,, 1L] - sigma[,, 2L])
  c(
    mean = n_cells * sum(d * solve(s_w, d)),
    cov = 0.5 * sum(diag(a %*% a))
  )
}

# The DWLS and lsei baselines read the bulk per cell, y / c.
lsei_plug_in <- function(y, mean_signature_matrix, mean_scale) {
  deconvolute_ratios_deconrnaseq(y / mean_scale, mean_signature_matrix)
}

#' DWLS on the bulk per cell
#'
#' One call of `DWLS::solveDampenedWLS()`, as `deconvolute_dwls()` in
#' DeCovarT_reproducibility (R/04_01_stats_deconvolution.R). Its dampening
#' constant is chosen on random halves of the genes, with the solver's own
#' fixed seeds; the caller's RNG state is restored and its console output
#' is dropped.
dwls_plug_in <- function(y, mean_signature_matrix, mean_scale) {
  bulk <- y / mean_scale
  utils::capture.output(
    estimate <- withr::with_preserve_seed(
      DWLS::solveDampenedWLS(S = mean_signature_matrix, B = bulk)
    )
  )
  stats::setNames(
    as.numeric(estimate[colnames(mean_signature_matrix)]),
    colnames(mean_signature_matrix)
  )
}

diagonal_sigma <- function(Sigma) {
  for (j in seq_len(dim(Sigma)[3L])) {
    Sigma[,, j] <- diag(diag(Sigma[,, j]))
  }
  Sigma
}

#' Methods and intervals of the bivariate single-cell benchmark
#'
#' @param intervals If `FALSE` (tests), no confidence intervals.
#' @param bootstrap If `TRUE`, Newton-Raphson only, with the bootstrap menu.
#' @param dwls If `FALSE`, DWLS is left out.
bivariate_single_cell_deconvolution_functions <- function(
  itmax = 100L,
  epsilon = 1e-6,
  intervals = TRUE,
  bootstrap = FALSE,
  n_boot = 199L,
  dwls = TRUE
) {
  solver <- list(epsilon = epsilon, itmax = itmax)
  menu <- function(methods, optimiser = "newton", ...) {
    if (!isTRUE(intervals)) {
      return(NULL)
    }
    list(methods = methods, optimiser = optimiser, ...)
  }
  arms <- c(oracle = "true", nominal = "nominal")
  if (isTRUE(bootstrap)) {
    newton <- function(scales, sdlog = 0) {
      list(
        FUN = deconvolute_single_cell,
        additional_parameters = c(solver, optimiser = "newton"),
        scales = scales,
        ci = list(
          methods = "bootstrap",
          optimiser = "newton",
          n_boot = n_boot,
          n_cells_sdlog = sdlog
        )
      )
    }
    return(list(
      "Newton-Raphson (oracle n)" = newton("true"),
      "Newton-Raphson (nominal n)" = newton("nominal"),
      # Relabelled `bootstrap_nprior` of the nominal arm after the run.
      "Newton-Raphson (nominal n, n prior)" = newton(
        "nominal",
        sqrt(log1p(FIG05_N_CELLS_CV^2))
      )
    ))
  }
  per_arm <- lapply(names(arms), function(arm) {
    scales <- arms[[arm]]
    decovart <- function(optimiser, ci_methods) {
      list(
        FUN = deconvolute_single_cell,
        additional_parameters = c(solver, optimiser = optimiser),
        scales = scales,
        ci = menu(ci_methods, optimiser)
      )
    }
    methods <- list(
      "lsei" = list(FUN = lsei_plug_in, scales = scales),
      "DWLS" = list(FUN = dwls_plug_in, scales = scales),
      "Newton-Raphson" = decovart(
        "newton",
        c("wald", "ilr", "profile", "oracle")
      ),
      "L-BFGS-B" = decovart("lbfgsb", c("wald", "ilr")),
      "Marquardt-Levenberg" = decovart("marquardt", c("wald", "ilr")),
      "Firth" = list(
        FUN = fit_firth_single_cell,
        scales = scales,
        ci = menu(c("wald", "ilr"))
      ),
      "DeCovarT diagonal" = c(
        decovart("newton", c("wald", "ilr")),
        list(working_sigma = diagonal_sigma)
      )
    )
    if (!isTRUE(dwls)) {
      methods$DWLS <- NULL
    }
    stats::setNames(methods, paste0(names(methods), " (", arm, " n)"))
  })
  c(
    unlist(per_arm, recursive = FALSE),
    list("NNLS" = list(FUN = deconvolute_ratios_nnls))
  )
}

#' Boundary likelihood-ratio test of p_2 = 0 for every saved bulk
#'
#' Newton-Raphson at known scales, oracle and nominal arms; chi-bar-square
#' calibration (half chi-square 0, half chi-square 1).
boundary_lrt_table <- function(benchmark, config, cores = 1L) {
  run_one <- function(i) {
    theta <- config$true_theta[[i]]
    data <- benchmark$data[[i]]
    nominal <- theta$n_cells
    purrr::map_dfr(seq_len(ncol(data$bulk)), function(s) {
      scales <- c(
        oracle = unname(data$nuisance["mean_scale", s]),
        nominal = nominal
      )
      purrr::map_dfr(names(scales), function(arm) {
        test <- tryCatch(
          suppressWarnings(lrt_single_cell(
            data$bulk[, s],
            theta$mu,
            theta$sigma,
            null_value = c(celltype_2 = 0),
            mean_scale = scales[[arm]],
            variance_scale = scales[[arm]]
          )),
          error = function(e) NULL
        )
        tibble::tibble(
          ID = config$ID[[i]],
          sample_id = s,
          arm = arm,
          p2_realised = data$truth[2L, s],
          statistic = if (is.null(test)) NA_real_ else test$statistic,
          p_value = if (is.null(test)) NA_real_ else test$p_value
        )
      })
    })
  }
  idx <- seq_len(nrow(config))
  out <- if (cores > 1L) {
    furrr::future_map(
      idx,
      run_one,
      .options = furrr::furrr_options(seed = TRUE)
    )
  } else {
    lapply(idx, run_one)
  }
  dplyr::bind_rows(out)
}

if (
  interactive() ||
    any(grepl(
      "fig05_bivariate_single_cell\\.R$",
      sub(
        "^--file=",
        "",
        grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
      )
    ))
) {
  # Run size and output. FIG05_DRY_RUN=1 runs N = 5 bulks into
  # output/fig05_dry/ to test the pipeline end to end;
  # FIG05_POSTPROCESS_ONLY=1 redraws the figures from the saved benchmark.
  DRY_RUN <- identical(Sys.getenv("FIG05_DRY_RUN"), "1")
  N_REPLICATES <- if (DRY_RUN) 5L else 200L
  N_REFERENCE <- if (DRY_RUN) 20L else 200L
  N_CELLS_REFERENCE <- 20L
  N_BOOT <- if (DRY_RUN) 19L else 199L
  CORES <- (parallel::detectCores() %/% 2L) %/% 3L
  POSTPROCESS_ONLY <- identical(Sys.getenv("FIG05_POSTPROCESS_ONLY"), "1")
  # DWLS chooses its dampening constant on random halves of the genes; with
  # G = 2 a half is one gene and DWLS 0.1.0 stops ("variable lengths
  # differ") on every bulk, so it is not run (TRUE to retry).
  RUN_DWLS <- FALSE
  OUT_DIR <- file.path("output", if (DRY_RUN) "fig05_dry" else "fig05")
  BENCHMARK_FILE <- file.path(OUT_DIR, "fig05_benchmark.rds")
  SEED <- 20261125L
  ITMAX <- 100L
  EPSILON <- 1e-6

  library(DeCovarT)
  local({
    a <- 0.5^abs(outer(seq_len(64L), seq_len(64L), "-"))
    if (
      max(abs(crossprod(chol(a)) - a)) > 1e-10 ||
        max(abs(solve(a) %*% a - diag(64L))) > 1e-8
    ) {
      stop(
        "BLAS returns wrong chol()/solve(); set OPENBLAS_CORETYPE=Haswell ",
        "before starting R.",
        call. = FALSE
      )
    }
  })
  source(file.path("scripts", "benchmark_metrics", "benchmark_artefacts.R"))
  source(
    file.path("scripts", "benchmark_metrics", "second_generation_methods.R")
  )
  source(file.path("scripts", "visualisation", "load_visualisation.R"))
  stopifnot(
    requireNamespace("MixSim", quietly = TRUE),
    !RUN_DWLS || requireNamespace("DWLS", quietly = TRUE)
  )

  cli::cli_h1("Figure 05: bivariate toy, single-cell model")

  if (!POSTPROCESS_ONLY) {
    set.seed(SEED)
    config <- build_bivariate_single_cell_config()
    cli::cli_alert_info(
      "{.val {nrow(config)}} scenarios, N = {.val {N_REPLICATES}}, {.val {CORES}} workers."
    )
    if (CORES > 1L) {
      future::plan(future::multisession, workers = CORES)
    }
    started <- proc.time()[["elapsed"]]
    out <- run_simulation_benchmark(
      scenario_config = config,
      deconvolution_functions = bivariate_single_cell_deconvolution_functions(
        itmax = ITMAX,
        epsilon = EPSILON,
        dwls = RUN_DWLS
      ),
      n = N_REPLICATES,
      cores = CORES,
      verbose = TRUE,
      likelihood = single_cell_likelihood(),
      return_data = TRUE,
      # Scenario i runs under seed SEED + i: the bootstrap pass below
      # redraws the same bulks.
      seed = SEED,
      reference = list(
        n_individuals = N_REFERENCE,
        n_cells = N_CELLS_REFERENCE,
        seed = SEED + 5000L
      )
    )

    cli::cli_alert_info("Parametric bootstrap on 8 scenarios.")
    boot_ids <- which(
      config$proportions %in%
        c("balanced", "moderately unbalanced") &
        config$n_cells %in% c(100, 1000) &
        paste(config$correlation_celltype1, config$correlation_celltype2) %in%
          c("0 0", "-0.8 0.8")
    )
    stopifnot(length(boot_ids) == 8L)
    boot <- run_simulation_benchmark(
      scenario_config = config,
      deconvolution_functions = bivariate_single_cell_deconvolution_functions(
        itmax = ITMAX,
        epsilon = EPSILON,
        bootstrap = TRUE,
        n_boot = N_BOOT
      ),
      n = N_REPLICATES,
      cores = CORES,
      verbose = TRUE,
      likelihood = single_cell_likelihood(),
      seed = SEED,
      scenarios = boot_ids
    )
    prior <- boot$intervals$algorithm == "Newton-Raphson (nominal n, n prior)"
    boot$intervals$ci_method[prior] <- "bootstrap_nprior"
    boot$intervals$algorithm[prior] <- "Newton-Raphson (nominal n)"
    out$intervals <- dplyr::bind_rows(out$intervals, boot$intervals)

    cli::cli_alert_info("RNA-Sieve and ReDeconv on the same bulks.")
    competitors <- run_second_generation_batch(
      competitor_scenarios(out, config),
      methods = c("RNA-Sieve", "ReDeconv"),
      workers = CORES
    )
    out <- score_competitors(
      out,
      competitors,
      config,
      likelihood = single_cell_likelihood()
    )

    cli::cli_alert_info("Boundary likelihood-ratio tests of p_2 = 0.")
    out$boundary_lrt <- boundary_lrt_table(out, config, cores = CORES)

    cli::cli_alert_success(
      "Benchmark finished in {.val {round(proc.time()[['elapsed']] - started)}} s."
    )
    save_benchmark(
      out,
      config,
      benchmark_meta(SEED, N_REPLICATES, N_REFERENCE),
      BENCHMARK_FILE
    )
  }

  benchmark <- readRDS(BENCHMARK_FILE)
  plot_data <- figure_plot_data(benchmark, figure = "fig05")
  saveRDS(plot_data, file.path(OUT_DIR, "fig05_plot_data.rds"))
  render_figure_books(plot_data, OUT_DIR)
  # Pages shown in vignettes/R3-bivariate-single-cell.qmd (pkgdown only).
  if (!DRY_RUN) {
    shown <- c(
      "density_visualisations/loglik_mean_shift_slice/balanced",
      "density_visualisations/loglik_correlation_slice/moderately_unbalanced",
      "performance_visualisations/heatmap_rmse/balanced_n_100",
      "performance_visualisations/heatmap_rmse/balanced_n_1000",
      "performance_visualisations/heatmap_rmse/moderately_unbalanced_n_1000",
      "performance_visualisations/ci_width_forest/balanced_n_1000_sb90_ba_n1e3",
      "performance_visualisations/ci_width_forest/balanced_n_100_sb9_ba_n1e2",
      "performance_visualisations/boundary_lrt"
    )
    target <- file.path("vignettes", "figures", "fig05")
    dir.create(target, showWarnings = FALSE, recursive = TRUE)
    file.copy(
      file.path(OUT_DIR, paste0(shown, ".svg")),
      file.path(target, paste0(gsub("/", "__", shown, fixed = TRUE), ".svg")),
      overwrite = TRUE
    )
  }
  cli::cli_alert_success("Figure 05 written to {.file {OUT_DIR}}.")
}
