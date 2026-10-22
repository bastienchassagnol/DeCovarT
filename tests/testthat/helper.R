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

# Three genes, three cell types with isotropic covariances: shared by the
# estimation (02_*) and inference (03_*) tests.
.inference_setup <- function() {
  signature <- matrix(
    c(20, 40, 15, 40, 20, 25, 25, 30, 35),
    nrow = 3,
    dimnames = list(paste0("gene_", 1:3), paste0("ct", 1:3))
  )
  Sigma <- array(
    0,
    dim = c(3, 3, 3),
    dimnames = list(
      rownames(signature),
      rownames(signature),
      colnames(signature)
    )
  )
  for (j in 1:3) {
    Sigma[,, j] <- diag(3) * c(1, 1.5, 1.2)[j]
  }
  list(signature = signature, Sigma = Sigma)
}

# Small single-cell bulk used by the interval tests (03_04).
.single_cell_interval_setup <- function() {
  mu <- matrix(
    c(2, 2, 0.2, 0.2, 0.2, 0.2, 3, 3),
    nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("ct", 1:2))
  )
  sigma <- array(
    0,
    dim = c(4, 4, 2),
    dimnames = list(rownames(mu), rownames(mu), colnames(mu))
  )
  sigma[,, 1] <- diag(0.05, 4)
  sigma[,, 2] <- diag(0.05, 4)
  y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
  list(mu = mu, sigma = sigma, y = y)
}

# Small simulation benchmarks shared by the visualisation tests (05_*).
.tiny_mc_benchmark <- function() {
  skip_on_os("windows")
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  theta <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
  scenario_config <- tibble::tibble(
    n_genes = c(2L, 2L),
    cosine = c(0, 0.5),
    true_theta = list(theta, theta)
  )
  withr::with_seed(
    11L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls)
      ),
      n = 6L,
      cores = 1L
    )
  )
}

.tiny_mc_benchmark_two_algos <- function() {
  skip_on_os("windows")
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  theta <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
  scenario_config <- tibble::tibble(
    n_genes = 2L,
    cosine = 0,
    true_theta = list(theta)
  )
  withr::with_seed(
    17L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls),
        "rlm" = list(FUN = deconvolute_ratios_rlm)
      ),
      n = 6L,
      cores = 1L
    )
  )
}
