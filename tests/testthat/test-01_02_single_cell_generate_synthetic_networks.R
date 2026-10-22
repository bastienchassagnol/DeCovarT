test_that("cellular_raw keeps sizes in the mean and the raw covariance", {
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
  expect_equal(ref$Sigma, sigma, ignore_attr = TRUE)
  expect_identical(ref$estimand, "cellular")
  expect_identical(ref$reference, "cellular_raw")
})

test_that("transcript_raw normalises the mean and divides the covariance", {
  mu <- cbind(c(2, 2, 0, 0), c(0, 0, 1, 1))
  sigma <- array(c(diag(4), 4 * diag(4)), dim = c(4, 4, 2))
  moments <- list(mean_profiles = mu, covariance_matrices = sigma)
  ref <- rescale_single_cell_moments(
    moments,
    cell_size = c(1, 5),
    reference = "transcript_raw"
  )

  expect_equal(colSums(ref$mean_signature_matrix), c(1, 1))
  expect_equal(ref$Sigma[,, 1], diag(4) / 1)
  expect_equal(ref$Sigma[,, 2], 4 * diag(4) / 5)
  expect_equal(colSums(ref$generative_mean), c(1, 5))
  expect_equal(ref$generative_sigma, sigma, ignore_attr = TRUE)
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
