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
      "intervals",
      "data",
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

test_that("the interval table reports every requested method", {
  genes <- paste0("g", 1:4)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 22, 15, 30, 22, 20, 25, 10),
    nrow = 4,
    dimnames = list(genes, cts)
  )
  sigma <- array(
    c(diag(4), 2 * diag(4)),
    dim = c(4, 4, 2),
    dimnames = list(genes, genes, cts)
  )
  config <- tibble::tibble(
    ID = "P1",
    true_theta = list(list(p = c(0.4, 0.6), mu = mu, sigma = sigma))
  )
  ci <- list(
    methods = c("wald", "ilr", "profile", "bootstrap", "oracle"),
    n_boot = 3L,
    n_draws = 100L
  )
  out <- withr::with_seed(
    1L,
    run_simulation_benchmark(
      config,
      list(
        "Newton-Raphson" = list(
          FUN = deconvolute_ratios_Newton_Raphson,
          ci = ci
        ),
        lsei = list(FUN = deconvolute_ratios_deconrnaseq)
      ),
      n = 3,
      return_data = TRUE,
      reference = list(n_individuals = 3L, seed = 10L)
    )
  )
  expect_setequal(out$intervals$ci_method, ci$methods)
  expect_true(all(out$intervals$algorithm == "Newton-Raphson"))
  expect_true(all(out$intervals$n_intervals == 3L))
  expect_true(all(out$intervals$mean_width > 0))
  expect_identical(dim(out$data[[1]]$bulk), c(4L, 3L))
  expect_identical(dim(out$data[[1]]$reference$reference$ct1), c(4L, 3L))
})

test_that("the single-cell model scores the realised composition and scales", {
  moments <- withr::with_seed(
    20261022L,
    simulate_hierarchical_grn_moments(
      n_genes = 10L,
      n_celltypes = 2L,
      mean_scale = 10,
      target_cosine = 0,
      nonnegative = TRUE,
      precision_shift = 0.5,
      precision_scale = 0.3,
      graph_model = "erdos_renyi"
    )
  )
  ref <- rescale_single_cell_moments(moments, cell_size = c(1, 5))
  config <- tibble::tibble(
    ID = "S1",
    true_theta = list(list(
      p = c(0.7, 0.3),
      mu = ref$mean_signature_matrix,
      sigma = ref$Sigma,
      n_cells = 1000,
      platform_scale = 5
    ))
  )
  out <- withr::with_seed(
    2L,
    run_simulation_benchmark(
      config,
      list(
        "Newton-Raphson" = list(
          FUN = deconvolute_single_cell,
          ci = list(methods = c("wald", "ilr_reml"), n_draws = 100L)
        )
      ),
      n = 3,
      likelihood = single_cell_likelihood(),
      return_data = TRUE
    )
  )
  expect_true(all(
    c(
      "mean_scale",
      "variance_scale",
      "platform_scale_a",
      "cell_count_n"
    ) %in%
      out$intervals$parameter
  ))
  expect_true(all(
    c("mean_scale", "variance_scale") %in% names(out$optimisation)
  ))
  truth <- out$data[[1]]$truth
  expect_false(isTRUE(all.equal(
    truth[, 1],
    c(0.7, 0.3),
    check.attributes = FALSE
  )))
  expect_equal(unname(out$data[[1]]$nuisance["mean_scale", ]), rep(5000, 3))
})

test_that("a scenario run alone with its seed reproduces the full run", {
  genes <- paste0("g", 1:3)
  cts <- c("ct1", "ct2")
  mu <- matrix(c(20, 22, 18, 21, 19, 23), 3, 2, dimnames = list(genes, cts))
  sigma <- array(
    c(diag(3), diag(3) * 2),
    dim = c(3, 3, 2),
    dimnames = list(genes, genes, cts)
  )
  config <- tibble::tibble(
    ID = c("A", "B"),
    true_theta = list(
      list(mu = mu, sigma = sigma, p = c(ct1 = 0.6, ct2 = 0.4)),
      list(mu = mu, sigma = sigma, p = c(ct1 = 0.3, ct2 = 0.7))
    )
  )
  fns <- list("NNLS" = list(FUN = deconvolute_ratios_nnls))
  # A simulator that fails on scenario B while the option is set.
  lik <- population_likelihood()
  simulate <- lik$simulate
  lik$simulate <- function(mean, Sigma, theta, n) {
    if (isTRUE(getOption("decovart.test_fail")) && theta$p[[1]] == 0.3) {
      stop("forced failure")
    }
    simulate(mean, Sigma, theta, n)
  }
  run <- function(...) {
    run_simulation_benchmark(config, fns, n = 4, likelihood = lik, seed = 7L, ...)
  }
  full <- run()

  withr::with_options(list(decovart.test_fail = TRUE), {
    expect_message(
      expect_error(run(), "forced failure"),
      "Scenario 2 failed \\(seed 9\\)"
    )
  })

  # Elapsed time and memory are measurements of the process, not results.
  results_only <- function(x) x[setdiff(names(x), c("elapsed_sec", "memory_bytes"))]
  alone <- run(scenarios = 2L)
  expect_equal(
    results_only(alone$optimisation),
    results_only(full$optimisation[full$optimisation$ID == "B", ]),
    ignore_attr = TRUE
  )
  expect_error(run(scenarios = 3L), "row indices")
})
