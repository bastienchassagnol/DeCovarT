spd_matrix <- function(G, seed) {
  withr::with_seed(seed, {
    a <- matrix(stats::rnorm(G * 2L * G), 2L * G, G)
    crossprod(a) / (2L * G) + diag(0.1, G)
  })
}

test_that("the Cholesky root reconstructs a definite and a semi-definite covariance", {
  sigma <- spd_matrix(6L, 1L)
  root <- .covariance_root(sigma)
  expect_equal(crossprod(root), sigma, tolerance = 1e-12)

  # Two identical genes: rank 2 of 3, chol() fails, the pivoted factor is
  # used.
  singular <- matrix(c(1, 1, 0, 1, 1, 0, 0, 0, 2), 3)
  expect_error(chol(singular))
  root_singular <- .covariance_root(singular)
  expect_equal(crossprod(root_singular), singular, tolerance = 1e-10)
})

test_that("Cholesky draws have the requested mean and covariance", {
  sigma <- spd_matrix(4L, 3L)
  mu <- c(a = 1, b = -2, c = 0.5, d = 3)
  draws <- withr::with_seed(4L, .rmvnorm_chol(20000L, mu, sigma))
  expect_identical(dim(draws), c(20000L, 4L))
  expect_identical(colnames(draws), names(mu))
  # Monte Carlo error of a mean is about sqrt(var / 20000) < 0.02 here.
  expect_lt(max(abs(colMeans(draws) - mu)), 0.05)
  expect_lt(max(abs(stats::cov(draws) - sigma)), 0.1 * max(diag(sigma)))
})

test_that("single-cell bulks at G = 200 are drawn without an eigendecomposition", {
  G <- 200L
  genes <- paste0("g", seq_len(G))
  cts <- c("ct1", "ct2")
  mu <- withr::with_seed(5L, matrix(
    stats::rexp(G * 2L, 1 / 5),
    G,
    2L,
    dimnames = list(genes, cts)
  ))
  sigma <- array(0, c(G, G, 2L), dimnames = list(genes, genes, cts))
  sigma[, , 1L] <- spd_matrix(G, 6L) * 1e-3
  sigma[, , 2L] <- spd_matrix(G, 7L) * 1e-3
  sim <- withr::with_seed(
    8L,
    simulate_bulk_from_cells(
      mu,
      sigma,
      p = c(0.7, 0.3),
      n_cells = 1000,
      n_samples = 200
    )
  )
  expect_identical(dim(sim$Y), c(G, 200L))
  expect_true(all(is.finite(sim$Y)))
  # E[y | counts] = mu %*% counts with a = 1.
  expected <- mu %*% sim$counts
  expect_lt(max(abs(rowMeans(sim$Y - expected))), 0.05 * max(mu))
})
