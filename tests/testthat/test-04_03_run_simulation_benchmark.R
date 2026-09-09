test_that("run_simulation_benchmark wraps simulate and deconvolute", {
  skip_on_os("windows")
  skip_if_not_installed("nnls")

  genes <- paste0("gene_", 1:2)
  cts <- paste0("celltype_", 1:2)
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(genes, cts)
  )
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  scenario_config <- tibble::tibble(
    ID = "B1_Ho",
    true_theta = list(list(
      p = c(0.5, 0.5),
      mu = mu,
      sigma = Sigma
    ))
  )

  out <- withr::with_seed(
    3L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls)
      ),
      n = 2L,
      cores = 1L
    )
  )

  expect_equal(nrow(out$config), 1L)
  expect_equal(nrow(out$optimisation), 2L)
  expect_true("elapsed_sec" %in% names(out$optimisation))
  expect_true("tv" %in% names(out$regression$global))
  expect_equal(out$config$nobservations[[1L]], 2L)
  expect_false(
    "parallel_scenarios" %in% names(formals(run_simulation_benchmark))
  )
  expect_named(
    out,
    c(
      "regression",
      "monte_carlo",
      "optimisation",
      "config",
      "theta_true",
      "descriptors",
      "supplementary",
      "call"
    )
  )
  expect_equal(nrow(out$descriptors), 1L)
  expect_true("f_cov" %in% names(out$descriptors))
  expect_true("mixsim_baromega" %in% names(out$descriptors))
  expect_true("hellinger" %in% names(out$descriptors))
  expect_equal(length(out$theta_true), 1L)
  expect_type(out$call, "language")
})

test_that("encode_bivariate_id records variance, entropy, and CLD", {
  expect_identical(
    DeCovarT:::.encode_bivariate_id(
      1L,
      "homoscedastic",
      "balanced",
      "small_CLD"
    ),
    "B1_Ho_Ba_Sm"
  )
  expect_identical(
    DeCovarT:::.encode_bivariate_id(
      12L,
      "heteroscedastic",
      "highly unbalanced",
      "large_CLD"
    ),
    "B12_He_Hi_Lg"
  )
})

test_that("slim_scenario_table drops aliases after canonicalising names", {
  tbl <- tibble::tibble(
    ID = "B1_Ho_Ba_Sm",
    scenario_idx = 1L,
    proportion_name = "balanced",
    centroid = "small_CLD",
    rho_ct1 = -0.8,
    correlation_celltype1 = -0.8
  )
  slim <- slim_scenario_table(tbl)
  expect_false("scenario_idx" %in% names(slim))
  expect_false("proportion_name" %in% names(slim))
  expect_false("centroid" %in% names(slim))
  expect_false("rho_ct1" %in% names(slim))
  expect_identical(slim$proportions, "balanced")
  expect_identical(slim$centroids, "small_CLD")
})

test_that("write_simulation_artefacts splits config, descriptors, and metrics", {
  skip_on_os("windows")
  skip_if_not_installed("nnls")

  genes <- paste0("gene_", 1:2)
  cts <- paste0("celltype_", 1:2)
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(genes, cts)
  )
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  scenario_config <- tibble::tibble(
    ID = "B1_Ho_Ba_Sm",
    scenario_idx = 1L,
    correlation_celltype1 = 0,
    correlation_celltype2 = 0,
    proportions = "balanced",
    variance = "homoscedastic",
    centroids = "small_CLD",
    true_theta = list(list(
      p = c(0.5, 0.5),
      mu = mu,
      sigma = Sigma
    ))
  )
  out <- withr::with_seed(
    3L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls)
      ),
      n = 2L,
      cores = 1L
    )
  )
  tmp <- withr::local_tempdir()
  paths <- write_simulation_artefacts(
    out,
    dir = tmp,
    stem = "bivariate",
    config = scenario_config
  )
  expect_true(file.exists(paths$config))
  expect_true(file.exists(paths$descriptors))
  expect_true(file.exists(paths$theta))
  expect_true(file.exists(paths$benchmark))
  cfg <- readRDS(paths$config)
  desc <- readRDS(paths$descriptors)
  metrics <- readRDS(paths$benchmark)
  expect_false("true_theta" %in% names(cfg))
  expect_false("scenario_idx" %in% names(cfg))
  expect_true("h_star" %in% names(desc))
  expect_true("mean_euclidean" %in% names(desc))
  expect_true("jeffreys" %in% names(desc))
  expect_false("config" %in% names(metrics))
  expect_true("optimisation" %in% names(metrics))
  assembled <- read_simulation_artefacts(tmp, "bivariate", assemble = TRUE)
  expect_equal(assembled$config$ID, "B1_Ho_Ba_Sm")
  expect_equal(nrow(assembled$optimisation), 2L)
})
