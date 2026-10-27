test_that("cellular_raw scales the mean by S_j and the covariance by S_j^2", {
  mu <- cbind(c(2, 2, 0, 0), c(0, 0, 1, 1))
  sigma <- array(c(diag(4), 4 * diag(4)), dim = c(4, 4, 2))
  moments <- list(mean_profiles = mu, covariance_matrices = sigma)
  ref <- rescale_single_cell_moments(
    moments,
    c(1, 5),
    reference = "cellular_raw"
  )

  expect_equal(colSums(ref$composition), c(1, 1))
  expect_equal(colSums(ref$generative_mean), c(1, 5))
  expect_equal(ref$mean_signature_matrix, ref$generative_mean)
  # Column totals were 4 and 2: Sigma_theta = Sigma / 16 and Sigma / 4.
  expect_equal(ref$generative_sigma[,, 1], diag(4) / 16, ignore_attr = TRUE)
  expect_equal(
    ref$generative_sigma[,, 2],
    25 * 4 * diag(4) / 4,
    ignore_attr = TRUE
  )
  expect_equal(ref$Sigma, ref$generative_sigma)
  expect_identical(ref$estimand, "cellular")
  expect_identical(ref$reference, "cellular_raw")
})

test_that("a common cell size leaves the standardised mean gap unchanged", {
  mu <- cbind(c(3, 1, 1, 1), c(1, 1, 3, 1))
  sigma <- array(rep(diag(4), 2), dim = c(4, 4, 2))
  moments <- list(mean_profiles = mu, covariance_matrices = sigma)
  mahalanobis_gap <- function(ref) {
    gap <- ref$generative_mean[, 1] - ref$generative_mean[, 2]
    drop(crossprod(gap, solve(ref$generative_sigma[,, 1], gap)))
  }
  expect_equal(
    mahalanobis_gap(rescale_single_cell_moments(moments, c(1, 1))),
    mahalanobis_gap(rescale_single_cell_moments(moments, c(5, 5)))
  )
})

test_that("transcript_raw keeps theta and S_j times the per-transcript covariance", {
  mu <- cbind(c(2, 2, 0, 0), c(0, 0, 1, 1))
  sigma <- array(c(diag(4), 4 * diag(4)), dim = c(4, 4, 2))
  moments <- list(mean_profiles = mu, covariance_matrices = sigma)
  ref <- rescale_single_cell_moments(
    moments,
    cell_size = c(1, 5),
    reference = "transcript_raw"
  )

  # Sigma_theta = Sigma / 16 and 4 Sigma / 4; the fit uses S_j Sigma_theta.
  expect_equal(colSums(ref$mean_signature_matrix), c(1, 1))
  expect_equal(ref$Sigma[,, 1], diag(4) / 16, ignore_attr = TRUE)
  expect_equal(ref$Sigma[,, 2], 5 * diag(4), ignore_attr = TRUE)
  expect_equal(colSums(ref$generative_mean), c(1, 5))
  expect_equal(ref$generative_sigma[,, 2], 25 * diag(4), ignore_attr = TRUE)
  expect_identical(ref$estimand, "transcript")
})

test_that("a scalar cell size is repeated across types", {
  mu <- cbind(c(1, 1, 1, 1), c(2, 0, 0, 0))
  sigma <- array(rep(diag(4), 2), dim = c(4, 4, 2))
  ref <- rescale_single_cell_moments(
    list(mean_profiles = mu, covariance_matrices = sigma),
    cell_size = 3
  )
  expect_equal(ref$cell_size, c(3, 3))
  expect_equal(colSums(ref$generative_mean), c(3, 3))
})

test_that("the scenario design is the proportion by cell-size grid", {
  design <- single_cell_scenario_design()
  expect_identical(design$n_genes, 10L)
  expect_identical(design$n_cells, 10000)
  expect_equal(design$platform_scale, c(1, 2))
  expect_named(
    design$scenarios,
    c(
      "balanced_equal_size",
      "balanced_varying_size",
      "imbalanced_equal_size",
      "imbalanced_varying_size"
    )
  )
  expect_equal(design$scenarios$imbalanced_equal_size$proportions, c(0.8, 0.2))
  expect_equal(design$scenarios$balanced_varying_size$cell_size, c(1, 5))
})

test_that("negative reference means are rejected", {
  mu <- cbind(c(1, -1, 0, 2), c(0, 1, 1, 0))
  sigma <- array(rep(diag(4), 2), dim = c(4, 4, 2))
  expect_error(
    rescale_single_cell_moments(
      list(mean_profiles = mu, covariance_matrices = sigma),
      c(1, 1)
    ),
    "nonnegative"
  )
})
