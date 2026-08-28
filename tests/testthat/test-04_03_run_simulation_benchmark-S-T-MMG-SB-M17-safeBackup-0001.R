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
  expect_equal(nrow(out$simulations), 2L)
  expect_true("model_mse" %in% names(out$simulations))
  expect_equal(out$config$nobservations[[1L]], 2L)
})

test_that(".default_parallel_cores uses half of detected cores", {
  half <- .default_parallel_cores()
  expect_type(half, "integer")
  expect_true(half >= 1L)
  expect_true(half <= max(1L, floor(parallel::detectCores() / 2L)))
})
