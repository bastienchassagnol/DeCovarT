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
  expect_equal(mc, mix, tolerance = 0.02)
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

test_that("AIRM distance is zero on the diagonal and inversion-invariant", {
  a <- matrix(c(2, 0.3, 0.3, 1.2), nrow = 2L)
  expect_equal(spd_affine_invariant_distance(a, a), 0, tolerance = 1e-10)
  b <- matrix(c(1.1, -0.2, -0.2, 2.4), nrow = 2L)
  expect_equal(
    spd_affine_invariant_distance(a, b),
    spd_affine_invariant_distance(solve(a), solve(b)),
    tolerance = 1e-8
  )
  theta <- list(
    mu = cbind(c(0, 0), c(1, 0)),
    sigma = array(c(a, b), dim = c(2, 2, 2))
  )
  expect_gt(compute_average_riemannian(theta), 0)
})

test_that("scaling covariances preserves precision zeros and orders overlap", {
  skip_if_not_installed("igraph")
  w <- matrix(
    c(
      0,
      1,
      0,
      0,
      1,
      0,
      1,
      0,
      0,
      1,
      0,
      1,
      0,
      0,
      1,
      0
    ),
    nrow = 4L
  )
  omega <- build_normalised_precision(w, precision_shift = 0.2)
  zero_mask <- abs(w) < .Machine$double.eps
  diag(zero_mask) <- FALSE
  expect_true(all(abs(omega[zero_mask]) < 1e-10))

  sigma_j <- solve(omega)
  sigma <- array(c(sigma_j, sigma_j), dim = c(4, 4, 2))
  mu <- cbind(c(4, 0, 0, 0), c(0, 4, 0, 0))
  p <- c(0.5, 0.5)
  targets <- c(0.05, 0.15, 0.30)
  calibrated <- lapply(targets, function(tgt) {
    scale_covariances_to_overlap(
      mu = mu,
      sigma = sigma,
      p = p,
      target = tgt,
      n_mc = 400L,
      s_range = c(1e-2, 80),
      seed = 11L,
      verbose = FALSE
    )
  })
  scales <- vapply(calibrated, `[[`, numeric(1), "scale")
  bars <- vapply(calibrated, `[[`, numeric(1), "baromega")
  expect_lt(bars[[1L]], bars[[3L]])
  expect_true(all(diff(scales) > -1e-8))
  omega_s <- solve(calibrated[[3L]]$sigma[,, 1L])
  expect_true(all(abs(omega_s[zero_mask]) < 1e-8))
})
