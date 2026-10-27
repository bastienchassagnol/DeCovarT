test_that("the reference panel has one column per individual and cell", {
  mu <- matrix(
    c(5, 1, 1, 1, 5, 1),
    nrow = 3,
    dimnames = list(paste0("g", 1:3), c("ct1", "ct2"))
  )
  sigma <- array(
    c(diag(3), 4 * diag(3)),
    dim = c(3, 3, 2),
    dimnames = list(rownames(mu), rownames(mu), colnames(mu))
  )
  panel <- simulate_reference_replicates(
    mu,
    sigma,
    n_individuals = 3,
    n_cells = 2,
    seed = 1
  )
  expect_named(panel$reference, c("ct1", "ct2"))
  expect_identical(dim(panel$reference$ct2), c(3L, 6L))
  expect_identical(panel$individual$ct1, rep(1:3, each = 2))
  expect_identical(rownames(panel$reference$ct1), rownames(mu))
})

test_that("the panel follows the cell-type moments and is reproducible", {
  mu <- matrix(
    c(5, 1, 1, 5),
    nrow = 2,
    dimnames = list(c("g1", "g2"), c("ct1", "ct2"))
  )
  sigma <- array(
    c(diag(2), 9 * diag(2)),
    dim = c(2, 2, 2),
    dimnames = list(rownames(mu), rownames(mu), colnames(mu))
  )
  panel <- simulate_reference_replicates(mu, sigma, 4000, seed = 2)
  expect_equal(rowMeans(panel$reference$ct1), mu[, 1], tolerance = 0.05)
  expect_equal(
    unname(diag(stats::cov(t(panel$reference$ct2)))),
    c(9, 9),
    tolerance = 0.1
  )
  expect_identical(
    simulate_reference_replicates(mu, sigma, 3, seed = 2)$reference,
    simulate_reference_replicates(mu, sigma, 3, seed = 2)$reference
  )
  expect_error(simulate_reference_replicates(mu, sigma, 1), "n_individuals")
})
