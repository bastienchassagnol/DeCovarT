test_that("scaling covariances preserves precision zeros and orders overlap", {
  skip_if_not_installed("igraph")
  skip_if_not(exists("scale_covariances_to_overlap", mode = "function"))
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
