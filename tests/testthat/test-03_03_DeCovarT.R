test_that("Use numDeriv to check numerically derivative values", {
  ## numDeriv::grad core options used below (method = "Richardson"):
  ## - eps: initial finite-difference step size (here 1e-4).
  ## - r: number of Richardson extrapolation iterations (here 6; default 4).
  ##   Larger r improves accuracy at the cost of more function evaluations
  ##   zero.tol, show.details. See ?numDeriv::grad.
  setup <- withr::with_seed(3L, {
    mean_signature_matrix <- matrix(c(20, 40, 40, 20), nrow = 2)
    p <- c(0.5, 0.5)
    num_genes <- nrow(mean_signature_matrix)
    num_celltypes <- ncol(mean_signature_matrix)
    y <- mean_signature_matrix %*% p + rnorm(nrow(mean_signature_matrix))
    Sigma <- array(
      c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
      dim = c(num_genes, num_genes, num_celltypes)
    )
    list(
      mean_signature_matrix = mean_signature_matrix,
      p = p,
      y = y,
      Sigma = Sigma
    )
  })

  ##----------------------------------------------------------------
  ##                  test gradient log-likelihood                 -
  ##----------------------------------------------------------------
  jacobian_mapping_numerical <- numDeriv::grad(
    loglik_multivariate,
    setup$p,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )

  jacobian_mapping_theoretical <- gradient_loglik_unconstrained(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    jacobian_mapping_numerical,
    jacobian_mapping_theoretical
  )
})


test_that("Vefiy that the DeCovarT algorithm is working, using a toy example", {
  # Define simulation parameters --------------------------------------------
  mean_signature_matrix <- matrix(
    c(20, 40, 40, 20),
    nrow = 2,
    dimnames = list(paste0("gene_", 1:2), paste0("celltype_", 1:2))
  )
  p <- c(0.5, 0.5)
  num_genes <- nrow(mean_signature_matrix)
  num_celltypes <- ncol(mean_signature_matrix)
  Sigma <- array(
    c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
    dim = c(num_genes, num_genes, num_celltypes),
    dimnames = list(
      paste0("gene_", 1:2),
      paste0("gene_", 1:2),
      paste0("celltype_", 1:2)
    )
  )

  # Simulate the bulk mixture as a convolution of bivariate Gaussian
  # distributions ----------------------------------------------------
  simulated_data <- withr::with_seed(
    3L,
    simulate_bulk_mixture(
      signature_matrix = mean_signature_matrix,
      Sigma,
      p = c(0.5, 0.5),
      n = 2000
    )
  )
  y_simu <- simulated_data$Y
  X_simu <- simulated_data$mean_signature_matrix

  # Deconvolution of synthetic data -------------------
  inferred_ratios <- deconvolute_ratios_Newton_Raphson(
    y = y_simu[, 2, drop = FALSE],
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    true_ratios = c(0.5, 0.5),
    epsilon = 10^-4,
    itmax = 200
  )
})









test_that("Benchmark standard deconvolution algorithms against DeCovarT", {
  deconvolution_functions <- list(
    "lm" = list(FUN = deconvolute_ratios_abbas),
    "nnls" = list(FUN = deconvolute_ratios_nnls),
    "lsei" = list(FUN = deconvolute_ratios_deconRNASeq),
    "LBFGS" = list(
      FUN = deconvolute_ratios_LBFGS,
      additional_parameters = list(epsilon = 10^-3, itmax = 200)
    ),
    "gradient" = list(
      FUN = deconvolute_ratios_first_order,
      additional_parameters = list(epsilon = 10^-3, itmax = 200)
    ),
    "Newton-Raphson" = list(
      FUN = deconvolute_ratios_second_order,
      additional_parameters = list(epsilon = 10^-3, itmax = 200)
    ),
    "Marquardt-Levenberg" = list(
      FUN = deconvolute_ratios_Marquardt_Levenberg,
      additional_parameters = list(epsilon = 10^-3, itmax = 200)
    ),
    "SA" = list(
      FUN = deconvolute_ratios_simulated_annealing,
      additional_parameters = list(epsilon = 10^-3, itmax = 200)
    )
  )

  bivariate_scenario <- withr::with_seed(
    seed = 3L,
    benchmark_bivariate_gaussian_convolutions(
      proportions = list("balanced" = c(0.50, 0.50)),
      n = 2,
      corr_sequence = c(-0.75, 0.75),
      signature_matrices = list(
        "small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)
      )
    ),
    .rng_kind = "L'Ecuyer-CMRG"
  )
  bivariate_configuraton <- readRDS(testthat::test_path(
    "fixtures",
    "bivariate_configuraton.rds"
  ))
  bivariate_estimation <- readRDS(testthat::test_path(
    "fixtures",
    "bivariate_estimation.rds"
  ))

  expect_equal(bivariate_configuraton, bivariate_scenario$config)
  expect_equal(bivariate_estimation, bivariate_scenario$simulations)
})
