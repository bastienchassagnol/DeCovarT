test_that("run_simulation_benchmark wraps simulate and deconvolute", {
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
