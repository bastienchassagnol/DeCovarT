# In-memory bivariate toy scenario for tests.
# See r-pkgs.org testing-advanced: create useful_things with a helper
# when construction is fiddly but cheap. Never write to the package
# tree, getwd(), or the user's home filespace (CRAN policy).

new_bivariate_toy_scenario <- function(
  seed = 3L,
  n = 2L,
  itmax = 10L,
  epsilon = 1e-3
) {
  pkg_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/"
  )
  scenario_script <- file.path(
    pkg_root,
    "scripts",
    "configure_bivariate_toy_scenarios.R"
  )
  if (!file.exists(scenario_script)) {
    return(NULL)
  }
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    return(NULL)
  }

  source(scenario_script, local = TRUE)

  scenario_config <- build_bivariate_scenario_config(
    proportions = list("balanced" = c(0.50, 0.50)),
    corr_sequence = 0,
    diagonal_terms = list("homoscedastic" = c(1, 1)),
    signature_matrices = list(
      "small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)
    )
  )
  deconvolution_functions <- bivariate_toy_deconvolution_functions(
    itmax = itmax,
    epsilon = epsilon
  )

  withr::with_seed(
    seed = seed,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = deconvolution_functions,
      n = n,
      cores = 1L
    )
  )
}
