###############################################################################
###  Rerun RNA-Sieve on the saved fig03 and fig06 benchmarks               ###
###############################################################################
#
# The fig03 (28 Oct 2026, 1c058d9) and fig06 (29 Oct 2026, f92538c) runs
# called the released rnasieve `_minimize_alpha`, which under cvxpy 1.6
# never moves the proportions from their start. The runner now installs the
# corrected update (scripts/benchmark_metrics/python/run_rnasieve.py). This
# script reruns RNA-Sieve only, on the saved bulks, replaces its rows in the
# benchmark (old file kept as `<file>.pre_rnasieve_fix.rds`) and redraws the
# figures. Other methods are untouched.
#
# The saved benchmarks drop the reference panels, so they are redrawn with
# the seeds of the original runs (SEED + 5000 + scenario index) and the
# sampler of those runs (MASS::mvrnorm; the package now uses a Cholesky
# draw, which gives other panels from the same seed).
#
#   export OPENBLAS_CORETYPE=Haswell OPENBLAS_NUM_THREADS=1
#   nohup Rscript scripts/rerun_rnasieve_fig03_fig06.R \
#     > "logs/rnasieve_rerun_$(date +%F).log" 2>&1 &

suppressPackageStartupMessages(library(DeCovarT))
source(file.path("scripts", "benchmark_metrics", "benchmark_artefacts.R"))
source(file.path("scripts", "benchmark_metrics", "second_generation_methods.R"))
source(file.path("scripts", "visualisation", "load_visualisation.R"))

CORES <- (parallel::detectCores() %/% 2L) %/% 3L
RUNS <- list(
  fig03 = list(
    seed = 20260807L,
    n_individuals = 200L,
    n_cells = 1L,
    likelihood = population_likelihood
  ),
  fig06 = list(
    seed = 20261024L,
    n_individuals = 200L,
    n_cells = 20L,
    likelihood = single_cell_likelihood
  )
)

# Reference panel exactly as simulate_reference_replicates() drew it at the
# commits of the original runs.
old_reference <- function(mu, Sigma, n_individuals, n_cells, seed) {
  n_columns <- n_individuals * n_cells
  individual <- rep(seq_len(n_individuals), each = n_cells)
  draw <- function(j) {
    x <- MASS::mvrnorm(n = n_columns, mu = mu[, j], Sigma = Sigma[,, j])
    x <- matrix(t(x), nrow = nrow(mu))
    dimnames(x) <- list(
      rownames(mu),
      paste0(
        "individual",
        individual,
        "_cell",
        sequence(rep(n_cells, n_individuals))
      )
    )
    x
  }
  reference <- withr::with_seed(seed, lapply(seq_len(ncol(mu)), draw))
  list(
    reference = stats::setNames(reference, colnames(mu)),
    individual = stats::setNames(
      rep(list(individual), ncol(mu)),
      colnames(mu)
    )
  )
}

drop_rnasieve <- function(tbl) tbl[tbl$algorithm != "RNA-Sieve", , drop = FALSE]

for (figure in names(RUNS)) {
  run <- RUNS[[figure]]
  out_dir <- file.path("output", figure)
  file <- file.path(out_dir, paste0(figure, "_benchmark.rds"))
  cli::cli_h1("{figure}: RNA-Sieve rerun")
  benchmark <- readRDS(file)
  config <- benchmark$config
  config$true_theta <- unname(benchmark$theta_true[config$ID])
  stopifnot(identical(names(benchmark$data), config$ID))

  data <- lapply(seq_len(nrow(config)), function(i) {
    theta <- config$true_theta[[i]]
    mu <- if (is.null(theta$generative_mean)) {
      theta$mu
    } else {
      theta$generative_mean
    }
    sigma <- if (is.null(theta$generative_sigma)) {
      theta$sigma
    } else {
      theta$generative_sigma
    }
    c(
      benchmark$data[[i]],
      list(
        reference = old_reference(
          mu,
          sigma,
          run$n_individuals,
          run$n_cells,
          run$seed + 5000L + i
        )
      )
    )
  })
  started <- proc.time()[["elapsed"]]
  competitors <- run_second_generation_batch(
    competitor_scenarios(list(data = data), config),
    methods = "RNA-Sieve",
    workers = CORES
  )
  cli::cli_alert_success(
    "RNA-Sieve done in {.val {round(proc.time()[['elapsed']] - started)}} s."
  )

  benchmark$optimisation <- drop_rnasieve(benchmark$optimisation)
  benchmark$intervals <- drop_rnasieve(benchmark$intervals)
  benchmark$monte_carlo <- drop_rnasieve(benchmark$monte_carlo)
  benchmark$regression$global <- drop_rnasieve(benchmark$regression$global)
  benchmark$regression$cell_type <- drop_rnasieve(
    benchmark$regression$cell_type
  )
  benchmark$data <- data
  benchmark <- score_competitors(
    benchmark,
    competitors,
    config,
    likelihood = run$likelihood()
  )

  file.copy(
    file,
    sub("\\.rds$", ".pre_rnasieve_fix.rds", file),
    overwrite = FALSE
  )
  benchmark$data <- lapply(benchmark$data, function(x) {
    x[c("bulk", "truth", "nuisance")]
  })
  names(benchmark$data) <- config$ID
  benchmark$meta$rnasieve_rerun <- list(
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    git_commit = tryCatch(
      system("git rev-parse --short HEAD", intern = TRUE),
      error = function(e) NA_character_
    )
  )
  saveRDS(benchmark, file)
  cli::cli_alert_success("Saved {.file {file}}.")

  benchmark <- readRDS(file)
  plot_data <- figure_plot_data(benchmark, figure = figure)
  saveRDS(plot_data, file.path(out_dir, paste0(figure, "_plot_data.rds")))
  render_figure_books(plot_data, out_dir)
  cli::cli_alert_success("{figure} figures redrawn.")
}
