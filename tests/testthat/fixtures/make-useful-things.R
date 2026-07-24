# Regenerates the bivariate fixtures used by
# tests/testthat/test-03_03_DeCovarT.R.
#
# Run from the package root after `devtools::load_all()`:
#   source("tests/testthat/fixtures/make-useful-things.R")

deconvolution_functions <- list(
  "nnls" = list(FUN = deconvolute_ratios_nnls),
  "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
  "LBFGS" = list(
    FUN = deconvolute_ratios_L_BFGS_B,
    additional_parameters = list(epsilon = 10^-3, itmax = 10)
  ),
  "gradient" = list(
    FUN = deconvolute_ratios_gradient_descent,
    additional_parameters = list(epsilon = 10^-3, itmax = 10)
  ),
  "Newton-Raphson" = list(
    FUN = deconvolute_ratios_Newton_Raphson,
    additional_parameters = list(epsilon = 10^-3, itmax = 10)
  ),
  "Marquardt-Levenberg" = list(
    FUN = deconvolute_ratios_Marquardt_Levenberg,
    additional_parameters = list(epsilon = 10^-3, itmax = 10)
  ),
  "SA" = list(
    FUN = deconvolute_ratios_simulated_annealing,
    additional_parameters = list(epsilon = 10^-3, itmax = 10)
  )
)

test_bivariate_scenario <- withr::with_seed(
  seed = 3L,
  benchmark_bivariate_gaussian_convolutions(
    proportions = list("balanced" = c(0.50, 0.50)),
    n = 2,
    corr_sequence = 0,
    diagonal_terms = list("homoscedastic" = c(1, 1)),
    signature_matrices = list(
      "small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)
    ),
    deconvolution_functions = deconvolution_functions,
    cores = 1
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
