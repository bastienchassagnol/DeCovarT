test_that("overlap_gaussian_mc matches MixSim BarOmega for J = 3, G = 2", {
  skip_if_not_installed("MixSim")
  theta <- list(
    p = rep(1 / 3, 3),
    mu = cbind(c(0, 0), c(3, 0), c(0, 3)),
    sigma = array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
  )
  mix <- MixSim::overlap(
    Pi = theta$p,
    Mu = t(theta$mu),
    S = theta$sigma
  )$BarOmega
  mc <- overlap_gaussian_mc(theta, n_mc = 8000L, seed = 1L)$BarOmega
  # Relative tolerance must exceed QMC Monte Carlo error when `qrng`
  # is unavailable (plain uniforms); 0.02 can fail on the boundary.
  expect_equal(mc, mix, tolerance = 0.03)
})

test_that("overlap_gaussian_mc is close to MixSim BarOmega for J = 3, G = 3", {
  skip_if_not_installed("MixSim")
  theta <- list(
    p = rep(1 / 3, 3),
    mu = cbind(c(0, 0, 0), c(3, 0, 0), c(0, 3, 0)),
    sigma = array(
      c(diag(3), diag(3), diag(3)),
      dim = c(3, 3, 3)
    )
  )
  mix <- MixSim::overlap(
    Pi = theta$p,
    Mu = t(theta$mu),
    S = theta$sigma
  )$BarOmega
  mc <- overlap_gaussian_mc(theta, n_mc = 6000L, seed = 2L)$BarOmega
  # BarOmega averages three pairwise overlaps; allow QMC error.
  expect_equal(mc, mix, tolerance = 0.06)
})

test_that("compute_average_overlap uses Monte Carlo when G >= 4", {
  theta <- list(
    p = c(0.5, 0.5),
    mu = cbind(rep(0, 4), c(2, 0, 0, 0)),
    sigma = array(c(diag(4), diag(4)), dim = c(4, 4, 2))
  )
  if (requireNamespace("MixSim", quietly = TRUE)) {
    testthat::local_mocked_bindings(
      overlap = function(...) {
        stop("MixSim::overlap should not run when G >= 4", call. = FALSE)
      },
      .package = "MixSim"
    )
  }
  ov <- withr::with_seed(
    2L,
    compute_average_overlap(
      theta,
      n_mc = 400L,
      seed = 2L,
      verbose = FALSE
    )
  )
  expect_type(ov, "double")
  expect_length(ov, 1L)
  expect_gte(ov, 0)
  expect_lte(ov, 1)
})

test_that("closer means raise Monte Carlo overlap when G >= 4", {
  far <- list(
    p = c(0.5, 0.5),
    mu = cbind(rep(0, 4), c(6, 0, 0, 0)),
    sigma = array(c(diag(4), diag(4)), dim = c(4, 4, 2))
  )
  near <- far
  near$mu[1L, 2L] <- 1
  ov_far <- withr::with_seed(
    3L,
    compute_average_overlap(
      far,
      n_mc = 600L,
      seed = 3L,
      verbose = FALSE
    )
  )
  ov_near <- withr::with_seed(
    3L,
    compute_average_overlap(
      near,
      n_mc = 600L,
      seed = 3L,
      verbose = FALSE
    )
  )
  expect_lt(ov_far, ov_near)
})

test_that("AIRM average is zero on identical covariances and inversion-invariant", {
  a <- matrix(c(2, 0.3, 0.3, 1.2), nrow = 2L)
  b <- matrix(c(1.1, -0.2, -0.2, 2.4), nrow = 2L)
  pair <- function(x, y) {
    list(mu = cbind(c(0, 0), c(1, 0)), sigma = array(c(x, y), dim = c(2, 2, 2)))
  }
  expect_equal(compute_average_riemannian(pair(a, a)), 0, tolerance = 1e-10)
  expect_equal(
    compute_average_riemannian(pair(a, b)),
    compute_average_riemannian(pair(solve(a), solve(b))),
    tolerance = 1e-8
  )
  expect_gt(compute_average_riemannian(pair(a, b)), 0)
})

test_that(".parse_true_theta accepts valid and sample-wise theta", {
  theta <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(3, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  parsed <- .parse_true_theta(theta, second_moment = "SIGMA")
  expect_identical(parsed$J, 2L)
  expect_identical(parsed$G, 2L)
  theta$p <- cbind(c(0.6, 0.4), c(0.4, 0.6))
  expect_equal(.parse_true_theta(theta)$p, c(0.5, 0.5))
  expect_error(
    .parse_true_theta(list(mu = matrix(1, 2, 2), p = c(0.5, 0.5))),
    "sigma"
  )
})

test_that("compute_average_overlap: closer means overlap more", {
  skip_if_not_installed("MixSim")
  set.seed(1)
  far <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(8, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  near <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(1, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )

  ov_far <- compute_average_overlap(far)
  ov_near <- compute_average_overlap(near)
  expect_type(ov_far, "double")
  expect_length(ov_far, 1L)
  expect_gte(ov_far, 0)
  expect_lte(ov_far, 1)
  expect_lt(ov_far, ov_near)
})

test_that("compute_average_overlap: J mismatch errors", {
  theta <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(1, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_error(compute_average_overlap(theta, J = 3L), "`J`")
})

test_that("compute_average_jeffreys: identical means give ~0", {
  theta <- list(
    mu = cbind(c(0, 0), c(0, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_equal(compute_average_jeffreys(theta), 0, tolerance = .tol_srr)
})

test_that("compute_average_jeffreys: mean shift increases divergence", {
  close <- list(
    mu = cbind(c(0, 0), c(0.5, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  far <- list(
    mu = cbind(c(0, 0), c(4, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_lt(compute_average_jeffreys(close), compute_average_jeffreys(far))
})

test_that("compute_average_jeffreys: equi-balanced default when p omitted", {
  theta <- list(
    mu = cbind(c(0, 0), c(2, 0), c(0, 2)),
    sigma = array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
  )
  j1 <- compute_average_jeffreys(theta)
  theta$p <- rep(1 / 3, 3)
  j2 <- compute_average_jeffreys(theta)
  expect_equal(j1, j2, tolerance = .tol_srr)
})
