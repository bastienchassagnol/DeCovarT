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
      z = isometric_log_ratio(p),
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
    setup$z,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  gradient_theoretical <- as.numeric(gradient_loglik_constrained(
    setup$z,
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
    setup$z,
    method = "Richardson",
    method.args = list(eps = 1e-4, r = 6),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix,
    Sigma = setup$Sigma
  )
  hessian_theoretical <- hessian_loglik_constrained(
    setup$z,
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


test_that("Helmert ILR maps are orthonormal and invert", {
  v <- helmert_basis(4L)
  expect_equal(crossprod(v), diag(3), tolerance = 1e-12)
  expect_equal(drop(crossprod(v, rep(1, 4))), rep(0, 3), tolerance = 1e-12)
  p <- c(0.2, 0.3, 0.1, 0.4)
  z <- isometric_log_ratio(p)
  expect_equal(isometric_logistic(z), p, tolerance = 1e-10)
  expect_equal(
    isometric_log_ratio(isometric_logistic(z)),
    z,
    tolerance = 1e-10
  )
})


test_that("Isometric logistic Jacobian and Hessian match numDeriv", {
  setup <- .decovart_deriv_setup()
  z <- setup$z

  jacobian_numerical <- numDeriv::jacobian(isometric_logistic, z)
  jacobian_theoretical <- jacobian_isometric_logistic(z)
  testthat::expect_equal(
    jacobian_numerical,
    jacobian_theoretical,
    tolerance = 1e-4
  )

  hessian_theoretical <- hessian_isometric_logistic(z)
  for (i in seq_along(setup$p)) {
    hessian_numerical <- numDeriv::hessian(
      function(zz) isometric_logistic(zz)[i],
      z
    )
    testthat::expect_equal(
      hessian_numerical,
      hessian_theoretical[,, i],
      tolerance = 1e-4
    )
  }
})


test_that("ILR and ALR delta-method covariances agree on the simplex", {
  setup <- .decovart_deriv_setup()
  v_ilr <- vcov_ilr_delta(
    setup$p,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  v_alr <- vcov_alr_delta(
    setup$p,
    setup$mean_signature_matrix,
    setup$Sigma
  )
  expect_equal(v_ilr, v_alr, tolerance = 1e-8)
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

test_that("starting_simplex supports barycentre, Dirichlet, and QP", {
  p_bar <- starting_simplex(3L)
  expect_equal(p_bar, rep(1 / 3, 3), tolerance = 1e-12)

  p_dir <- withr::with_seed(11L, starting_simplex(3L, "dirichlet"))
  expect_length(p_dir, 3L)
  expect_equal(sum(p_dir), 1, tolerance = 1e-8)
  expect_true(all(p_dir > 0))

  skip_if_not_installed("limSolve")
  setup <- .decovart_deriv_setup()
  p_qp <- starting_simplex(
    ncol(setup$mean_signature_matrix),
    "qp",
    colnames(setup$mean_signature_matrix),
    y = setup$y,
    mean_signature_matrix = setup$mean_signature_matrix
  )
  expect_equal(sum(p_qp), 1, tolerance = 1e-8)
  expect_true(all(p_qp > 0))
})

test_that("Benchmark standard deconvolution algorithms against DeCovarT", {
  skip_on_os("windows")
  skip_if_not_installed("nnls")
  skip_if_not_installed("limSolve")
  skip_if_not_installed("MixSim")

  bivariate_scenario <- new_bivariate_toy_scenario()
  skip_if(
    is.null(bivariate_scenario),
    "scripts/ is not in the package tarball"
  )

  expect_type(bivariate_scenario, "list")
  expect_named(
    bivariate_scenario,
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
  expect_equal(nrow(bivariate_scenario$config), 1L)
  expect_true(nrow(bivariate_scenario$optimisation) > 0L)
  expect_true("tv" %in% names(bivariate_scenario$regression$global))
})
