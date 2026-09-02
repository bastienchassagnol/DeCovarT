# Reproduce the bivariate toy benchmark (reduced grid for smoke testing).
#
# Full factorial grids and figure export live in scripts/bivariate_toy_model.R.
# Scenario parameters are built in scripts/configure_bivariate_toy_scenarios.R.
#
# Run from the package root:
#   Rscript scripts/run_bivariate_toy_benchmark.R

if (!requireNamespace("DeCovarT", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("Install DeCovarT or devtools before running this script.")
  }
}

source("scripts/configure_bivariate_toy_scenarios.R", local = TRUE)

scenario_config <- build_bivariate_scenario_config(
  proportions = list("balanced" = c(0.5, 0.5)),
  signature_matrices = list(
    "small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)
  ),
  corr_sequence = 0,
  diagonal_terms = list("homoscedastic" = c(1, 1))
)

deconvolution_functions <- bivariate_toy_deconvolution_functions(
  itmax = 10L,
  epsilon = 1e-3
)

set.seed(3L)
benchmark_out <- run_simulation_benchmark(
  scenario_config = scenario_config,
  deconvolution_functions = deconvolution_functions,
  n = 2L,
  cores = 1L
)

message(
  "Bivariate toy benchmark: ",
  nrow(benchmark_out$config),
  " scenario(s), ",
  nrow(benchmark_out$optimisation),
  " estimation row(s)."
)

out_dir <- file.path("simulations", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  benchmark_out$optimisation,
  file.path(out_dir, "bivariate_scenario.rds")
)
saveRDS(
  benchmark_out$config,
  file.path(out_dir, "complete_bivariate_configuration.rds")
)

message("Saved RDS artefacts under ", normalizePath(out_dir))
