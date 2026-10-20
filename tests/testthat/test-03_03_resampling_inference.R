# Resampling inference: parametric bootstrap, reference-sample bootstrap
# and the label-equivariance check.
#
#' @srrstats {RE4.3} Bootstrap calibrations (`bootstrap_decovart()`,
#'   `reference_bootstrap_decovart()`) are exercised below.

test_that("bootstrap_decovart reproduces the atom at zero under H0", {
  setup <- .inference_setup()
  replicates <- withr::with_seed(
    11L,
    simulate_bulk_mixture(
      signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      p = c(0.5, 0.5, 0),
      n = 4L
    )$Y
  )

  calibration <- withr::with_seed(
    5L,
    bootstrap_decovart(
      y = replicates,
      mean_signature_matrix = setup$signature,
      Sigma = setup$Sigma,
      null_value = c(ct3 = 0),
      n_boot = 60L
    )
  )

  expect_length(calibration$null_statistics, 60L)
  expect_true(all(calibration$null_statistics >= 0))
  # Chernoff / Self-Liang: about half the replicates land on the face.
  expect_gt(mean(calibration$null_statistics < 1e-8), 0.3)
  expect_lt(mean(calibration$null_statistics < 1e-8), 0.7)
  expect_gt(calibration$p_value, 0.05)
})

test_that("equivariance_check_decovart recovers the relabelled MLE", {
  setup <- .inference_setup()
  y <- drop(setup$signature %*% c(0.7, 0.2, 0.1))
  perm <- c(3L, 1L, 2L)

  check <- equivariance_check_decovart(
    y,
    setup$signature,
    setup$Sigma,
    perm = perm,
    itmax = 80
  )
  expect_equal(check$p_star, check$p_expected, tolerance = 1e-5)
  expect_lt(abs(check$loglik_diff), 1e-6)
  expect_identical(check$perm, perm)
})

test_that("reference_bootstrap_decovart resamples cells within type", {
  setup <- .inference_setup()
  p <- c(0.55, 0.30, 0.15)
  refs <- lapply(seq_len(3), function(j) {
    draws <- MASS::mvrnorm(
      n = 10L,
      mu = setup$signature[, j],
      Sigma = setup$Sigma[,, j]
    )
    out <- t(draws)
    rownames(out) <- rownames(setup$signature)
    out
  })
  names(refs) <- colnames(setup$signature)
  y <- drop(setup$signature %*% p)

  expect_error(
    reference_bootstrap_decovart(y, refs, n_boot = 2L, itmax = 40),
    "donor_ids"
  )

  boot <- withr::with_seed(
    13L,
    reference_bootstrap_decovart(
      y,
      refs,
      method = "cells",
      n_boot = 12L,
      itmax = 80
    )
  )
  expect_identical(boot$method, "cells")
  expect_identical(dim(boot$estimates), c(3L, 12L))
  expect_named(boot$p_hat, colnames(setup$signature))
  expect_equal(colSums(boot$estimates), rep(1, 12L), tolerance = 1e-6)
  expect_identical(dim(boot$interval), c(3L, 2L))
})

test_that("reference_bootstrap_decovart resamples donors within type", {
  setup <- .inference_setup()
  p <- c(0.55, 0.30, 0.15)
  refs <- lapply(seq_len(3), function(j) {
    draws <- MASS::mvrnorm(
      n = 12L,
      mu = setup$signature[, j],
      Sigma = setup$Sigma[,, j]
    )
    out <- t(draws)
    rownames(out) <- rownames(setup$signature)
    out
  })
  names(refs) <- colnames(setup$signature)
  donor_ids <- lapply(refs, function(x) {
    rep(paste0("d", 1:3), length.out = ncol(x))
  })
  y <- drop(setup$signature %*% p)

  boot <- withr::with_seed(
    14L,
    reference_bootstrap_decovart(
      y,
      refs,
      donor_ids = donor_ids,
      method = "donors",
      n_boot = 8L,
      itmax = 80
    )
  )
  expect_identical(boot$method, "donors")
  expect_identical(dim(boot$estimates), c(3L, 8L))
  expect_equal(colSums(boot$estimates), rep(1, 8L), tolerance = 1e-6)
})
