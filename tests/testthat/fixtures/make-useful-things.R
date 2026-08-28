# Regenerates the bivariate fixtures used by
# tests/testthat/test-03_03_DeCovarT.R.
#
# Run from the package root after `devtools::load_all()`:
#   source("tests/testthat/fixtures/make-useful-things.R")

source(
  file.path("scripts", "configure_bivariate_toy_scenarios.R"),
  local = TRUE
)

scenario_config <- build_bivariate_scenario_config(
  proportions = list("balanced" = c(0.50, 0.50)),
  corr_sequence = 0,
  diagonal_terms = list("homoscedastic" = c(1, 1)),
  signature_matrices = list(
    "small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)
  )
)

deconvolution_functions <- bivariate_toy_deconvolution_functions(
  itmax = 10L,
  epsilon = 1e-3
)

test_bivariate_scenario <- withr::with_seed(
  seed = 3L,
  run_simulation_benchmark(
    scenario_config = scenario_config,
    deconvolution_functions = deconvolution_functions,
    n = 2L,
    cores = 1L,
    parallel_scenarios = FALSE
  )
)

fixtures_dir <- if (requireNamespace("testthat", quietly = TRUE)) {
  testthat::test_path("fixtures")
} else {
  "tests/testthat/fixtures"
}
dir.create(fixtures_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  test_bivariate_scenario$config,
  file.path(fixtures_dir, "bivariate_configuration.rds")
)
saveRDS(
  test_bivariate_scenario$simulations,
  file.path(fixtures_dir, "bivariate_estimation.rds")
)

message("Wrote fixtures to ", normalizePath(fixtures_dir))
