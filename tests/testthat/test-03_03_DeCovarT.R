# Shared scenario (3 genes, 3 cell types) for the derivative checks below.
# A three-component simplex yields a non-trivial (2 x 2) constrained Hessian.
.decovart_deriv_setup <- function() {
  withr::with_seed(7L, {
    mean_signature_matrix <- matrix(
      c(20, 40, 15, 40, 20, 25, 25, 30, 35),
      nrow = 3,
      dimnames = list(paste0("gene_", 1:3), paste0("celltype_", 1:3))
    )
    Sigma <- array(0, dim = c(3, 3, 3))
    Sigma[,, 1] <- matrix(c(1, 0.3, 0.1, 0.3, 1, 0.2, 0.1, 0.2, 1), 3, 3)
    Sigma[,, 2] <- matrix(c(2, -0.2, 0, -0.2, 1.5, 0.3, 0, 0.3, 1), 3, 3)
    Sigma[,, 3] <- matrix(
      c(1.2, 0.1, 0.4, 0.1, 0.8, -0.1, 0.4, -0.1, 1.1),
      3,
      3
    )
    p <- c(0.2, 0.3, 0.5)
    y <- mean_signature_matrix %*% p + rnorm(3)
    list(
      mean_signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      p = p,
      rho = additive_log_ratio(p),
      y = y
    )
  })
}


test_that(".inner_product matches the closed-form quadratic product", {
  x <- c(1, 2)
  y <- c(3, -1)
  a <- matrix(c(2, 0.5, 0.5, 3), 2L)

  expect_equal(
    .inner_product(x, a, y),
    as.numeric(t(x) %*% a %*% y),
    tolerance = .tol_srr
  )
  expect_equal(
    .inner_product(x, a),
    as.numeric(t(x) %*% a %*% x),
    tolerance = .tol_srr
  )
})

test_that(".sigma_p_factorisation caches by exact (p, Sigma)", {
  setup <- .decovart_deriv_setup()
  first <- .sigma_p_factorisation(setup$p, setup$Sigma)
  second <- .sigma_p_factorisation(setup$p, setup$Sigma)

  expect_identical(first, second)
  expect_named(first, c("matrix", "chol", "log_det", "inverse"))
  expect_equal(
    first$matrix,
    .compute_global_variance(setup$p, setup$Sigma),
    tolerance = 1e-10
  )
  expect_equal(
    first$inverse,
    solve(first$matrix),
    tolerance = 1e-8
  )
})


test_that("Unconstrained objective gradient/Hessian match numDeriv", {
  setup <- .decovart_deriv_setup()

  gradient_numerical <- numDeriv::grad(
    loglik_multivariate,
    setup$p,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  gradient_theoretical <- gradient_loglik_unconstrained(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    gradient_numerical,
    gradient_theoretical,
    tolerance = 1e-4
  )

  hessian_numerical <- numDeriv::hessian(
    loglik_multivariate,
    setup$p,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  hessian_theoretical <- hessian_loglik_unconstrained(
    setup$p,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    hessian_numerical,
    hessian_theoretical,
    tolerance = 1e-4
  )
})


test_that("Constrained gradient and Hessian of the objective match numDeriv", {
  setup <- .decovart_deriv_setup()

  gradient_numerical <- numDeriv::grad(
    loglik_multivariate_constrained,
    setup$rho,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  gradient_theoretical <- as.numeric(gradient_loglik_constrained(
    setup$rho,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  ))
  testthat::expect_equal(
    gradient_numerical,
    gradient_theoretical,
    tolerance = 1e-4
  )

  hessian_numerical <- numDeriv::hessian(
    loglik_multivariate_constrained,
    setup$rho,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  hessian_theoretical <- hessian_loglik_constrained(
    setup$rho,
    setup$y,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  testthat::expect_equal(
    hessian_numerical,
    hessian_theoretical,
    tolerance = 1e-4
  )
})


test_that("Additive logistic Jacobian and Hessian match numDeriv", {
  setup <- .decovart_deriv_setup()

  jacobian_numerical <- numDeriv::jacobian(additive_logistic, setup$rho)
  jacobian_theoretical <- jacobian_additive_logistic(setup$rho)
  testthat::expect_equal(
    jacobian_numerical,
    jacobian_theoretical,
    tolerance = 1e-4
  )

  # Hessian tensor is stored slice-wise as [ , , i] = d^2 p_i / (d rho d rho^T)
  hessian_theoretical <- hessian_additive_logistic(setup$rho)
  for (i in seq_along(setup$p)) {
    hessian_numerical <- numDeriv::hessian(
      function(rho) additive_logistic(rho)[i],
      setup$rho
    )
    testthat::expect_equal(
      hessian_numerical,
      hessian_theoretical[,, i],
      tolerance = 1e-4
    )
  }
})


test_that("All deconvolute_ratios_* solvers return a valid simplex", {
  setup <- .decovart_deriv_setup()
  expect_valid_simplex <- function(estimated_ratios, label) {
    testthat::expect_length(estimated_ratios, length(setup$p))
    testthat::expect_false(anyNA(estimated_ratios), info = label)
    testthat::expect_true(
      all(estimated_ratios >= 0 & estimated_ratios <= 1),
      info = label
    )
    testthat::expect_equal(
      sum(estimated_ratios),
      1,
      tolerance = 1e-6,
      info = label
    )
  }

  # ------------------------------------------------------------------ #
  # Marquardt–Levenberg (damped Newton on unconstrained rho)           #
  # ------------------------------------------------------------------ #
  inferred_ml <- withr::with_seed(
    3L,
    suppressWarnings(
      deconvolute_ratios_Marquardt_Levenberg(
        y = setup$y,
        mean_signature_matrix = setup$mean_signature_matrix,
        Sigma = setup$Sigma,
        epsilon = 10^-4,
        itmax = 200
      )
    )
  )
  expect_valid_simplex(as.numeric(inferred_ml), "ML")

  # ------------------------------------------------------------------ #
  # Simulated annealing (zeroth-order on unconstrained rho)            #
  # ------------------------------------------------------------------ #
  inferred_sa <- withr::with_seed(
    3L,
    suppressWarnings(
      deconvolute_ratios_simulated_annealing(
        y = setup$y,
        mean_signature_matrix = setup$mean_signature_matrix,
        Sigma = setup$Sigma,
        epsilon = 10^-4,
        itmax = 200
      )
    )
  )
  expect_valid_simplex(as.numeric(inferred_sa), "SANN")

  # ------------------------------------------------------------------ #
  # L-BFGS-B (first-order, box-constrained in p-space)                 #
  # ------------------------------------------------------------------ #
  inferred_lbfgs <- withr::with_seed(
    3L,
    suppressWarnings(
      deconvolute_ratios_L_BFGS_B(
        y = setup$y,
        mean_signature_matrix = setup$mean_signature_matrix,
        Sigma = setup$Sigma,
        epsilon = 10^-4,
        itmax = 200
      )
    )
  )
  expect_valid_simplex(as.numeric(inferred_lbfgs), "L-BFGS-B")

  # ------------------------------------------------------------------ #
  # Newton–Raphson (second-order on unconstrained rho)                 #
  # ------------------------------------------------------------------ #
  inferred_nr <- withr::with_seed(
    3L,
    suppressWarnings(
      deconvolute_ratios_Newton_Raphson(
        y = setup$y,
        mean_signature_matrix = setup$mean_signature_matrix,
        Sigma = setup$Sigma,
        epsilon = 10^-4,
        itmax = 200
      )
    )
  )
  expect_valid_simplex(as.numeric(inferred_nr), "Newton")

  # ------------------------------------------------------------------ #
  # BFGS / gradient ascent (first-order on unconstrained rho)          #
  # ------------------------------------------------------------------ #
  inferred_gd <- withr::with_seed(
    3L,
    suppressWarnings(
      deconvolute_ratios_gradient_descent(
        y = setup$y,
        mean_signature_matrix = setup$mean_signature_matrix,
        Sigma = setup$Sigma,
        epsilon = 10^-4,
        itmax = 200
      )
    )
  )
  expect_valid_simplex(as.numeric(inferred_gd), "BFGS")
})

test_that("Benchmark standard deconvolution algorithms against DeCovarT", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("limSolve")
  pkg_root <- normalizePath(
    file.path(testthat::test_path(), "..", ".."),
    winslash = "/"
  )
  scenario_script <- file.path(
    pkg_root,
    "scripts",
    "configure_bivariate_toy_scenarios.R"
  )
  skip_if_not(
    file.exists(scenario_script),
    "scripts/ is not in the package tarball"
  )
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
    itmax = 10L,
    epsilon = 1e-3
  )

  bivariate_scenario <- withr::with_seed(
    seed = 3L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = deconvolution_functions,
      n = 2L,
      cores = 1L,
      parallel_scenarios = FALSE
    )
  )

  bivariate_configuration <- readRDS(testthat::test_path(
    "fixtures",
    "bivariate_configuration.rds"
  ))
  bivariate_estimation <- readRDS(testthat::test_path(
    "fixtures",
    "bivariate_estimation.rds"
  ))

  expect_equal(bivariate_configuration, bivariate_scenario$config)
  expect_equal(bivariate_estimation, bivariate_scenario$simulations)
})
