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
      z = isometric_log_ratio(p),
      y = y
    )
  })
}


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


test_that("Unconstrained Hessian is symmetric and matches numDeriv for J = 4", {
  # Regression: the former symmetrisation H[lower.tri] <- H[upper.tri]
  # paired the wrong entries as soon as J >= 4.
  withr::local_seed(3L)
  n_genes <- 8L
  n_celltypes <- 4L
  mu <- matrix(stats::rexp(n_genes * n_celltypes, 0.1), n_genes, n_celltypes)
  Sigma <- array(0, dim = c(n_genes, n_genes, n_celltypes))
  for (j in seq_len(n_celltypes)) {
    a <- matrix(stats::rnorm(n_genes^2), n_genes)
    Sigma[,, j] <- crossprod(a) / n_genes + diag(0.3, n_genes)
  }
  p <- c(0.4, 0.3, 0.2, 0.1)
  y <- drop(mu %*% p) + stats::rnorm(n_genes)

  hessian_theoretical <- hessian_loglik_unconstrained(p, y, mu, Sigma)
  hessian_numerical <- numDeriv::hessian(
    loglik_multivariate,
    p,
    y = y,
    mean_signature_matrix = mu,
    Sigma = Sigma
  )
  expect_true(isSymmetric(hessian_theoretical))
  expect_equal(hessian_theoretical, hessian_numerical, tolerance = 1e-6)
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


test_that("loglik_multivariate is the Gaussian log-density up to a constant", {
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  y <- drop(setup$signature %*% p) + c(0.4, -0.2, 0.1)

  covariance <- .compute_global_variance(p, setup$Sigma)
  residual <- y - drop(setup$signature %*% p)
  reference <- -0.5 *
    as.numeric(determinant(covariance, logarithm = TRUE)$modulus) -
    0.5 * drop(t(residual) %*% solve(covariance) %*% residual)

  expect_equal(
    loglik_multivariate(p, y, setup$signature, setup$Sigma),
    reference,
    tolerance = 1e-10
  )
})

test_that("loglik_multivariate matches dmvnorm up to the 2 pi constant", {
  skip_if_not_installed("mvtnorm")
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  y <- drop(setup$signature %*% p) + c(0.4, -0.2, 0.1)
  covariance <- .compute_global_variance(p, setup$Sigma)
  n_genes <- length(y)

  expect_equal(
    loglik_multivariate(p, y, setup$signature, setup$Sigma) -
      0.5 * n_genes * log(2 * pi),
    as.numeric(mvtnorm::dmvnorm(
      y,
      mean = drop(setup$signature %*% p),
      sigma = covariance,
      log = TRUE
    )),
    tolerance = 1e-10
  )
})

test_that("Cholesky backsolve matches the explicit-inverse Mahalanobis term", {
  setup <- .inference_setup()
  p <- c(0.45, 0.35, 0.20)
  y <- drop(setup$signature %*% p) + c(0.1, -0.3, 0.2)
  residual <- as.numeric(y - drop(setup$signature %*% p))
  sigma_p <- .sigma_p_factorisation(p, setup$Sigma)
  z <- backsolve(sigma_p$chol, residual, transpose = TRUE)

  expect_equal(
    sum(z * z),
    drop(crossprod(residual, sigma_p$inverse %*% residual)),
    tolerance = 1e-10
  )
})

test_that("cell-type permutation leaves the Gaussian log-likelihood invariant", {
  setup <- .inference_setup()
  p <- c(0.5, 0.3, 0.2)
  perm <- c(3L, 1L, 2L)
  y <- drop(setup$signature %*% p)
  mu_perm <- setup$signature[, perm, drop = FALSE]
  colnames(mu_perm) <- colnames(setup$signature)
  sigma_perm <- setup$Sigma[,, perm, drop = FALSE]
  if (!is.null(dimnames(setup$Sigma))) {
    dimnames(sigma_perm) <- dimnames(setup$Sigma)
  }

  expect_equal(
    loglik_multivariate(p, y, setup$signature, setup$Sigma),
    loglik_multivariate(p[perm], y, mu_perm, sigma_perm),
    tolerance = 1e-10
  )

  gene_perm <- c(3L, 1L, 2L)
  y_gene <- y[gene_perm]
  expect_false(
    isTRUE(all.equal(
      loglik_multivariate(p, y, setup$signature, setup$Sigma),
      loglik_multivariate(p, y_gene, setup$signature, setup$Sigma),
      tolerance = 1e-8
    ))
  )
})

test_that("isometric_logistic stays finite for extreme ILR coordinates", {
  p <- isometric_logistic(c(800, 1200))

  expect_false(anyNA(p))
  expect_equal(sum(p), 1, tolerance = 1e-12)
})

#' @srrstats {G5.9b} `multistart_decovart()` sequences several
#'   Dirichlet starts; spread of attained log-likelihoods is recorded.
test_that("multistart_decovart explores several starts", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.45, 0.35, 0.20))

  restarts <- withr::with_seed(
    3L,
    multistart_decovart(
      y = y,
      mean_signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      n_starts = 3L,
      itmax = 100
    )
  )

  expect_identical(dim(restarts$starts), c(3L, 4L))
  expect_length(restarts$logliks, 4L)
  expect_type(restarts$multimodal, "logical")
  expect_equal(sum(restarts$coefficients), 1, tolerance = 1e-6)
})


test_that("repair_simplex renormalises and rejects invalid inputs", {
  p_near <- c(0.2, 0.3, 0.5 + 1e-12)
  expect_equal(
    repair_simplex(p_near),
    p_near / sum(p_near),
    tolerance = .tol_srr
  )
  expect_equal(sum(repair_simplex(c(2, 2, 0))), 1, tolerance = .tol_srr)
  expect_error(repair_simplex(numeric()), "`p`")
  expect_error(repair_simplex(c(0.5, NA)), "`p`")
  expect_error(repair_simplex(c(-1, 2)), "negative")
  expect_error(repair_simplex(c(0, 0)), "positive")
})
