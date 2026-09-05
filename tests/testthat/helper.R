# In-memory bivariate toy scenario for tests.
# See r-pkgs.org testing-advanced: create useful_things with a helper
# when construction is fiddly but cheap. Never write to the package
# tree, getwd(), or the user's home filespace (CRAN policy).
#
# The fig02 script is copied into a withr temp directory and sourced
# from there so builders load without running the 972-scenario pipeline
# (that path is gated on Rscript --file=...fig02_bivariate_toy.R or an
# interactive session). The temp copy is deleted when this helper
# returns.

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
    "fig02_bivariate_toy.R"
  )
  if (!file.exists(scenario_script)) {
    return(NULL)
  }
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    return(NULL)
  }

  tmp_dir <- withr::local_tempdir()
  tmp_script <- file.path(tmp_dir, "fig02_bivariate_toy.R")
  file.copy(scenario_script, tmp_script)
  source(tmp_script, local = TRUE)

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
