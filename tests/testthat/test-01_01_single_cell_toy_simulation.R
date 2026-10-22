test_that("multinomial cell counts match the returned technical factors", {
  genes <- paste0("g", 1:4)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(2, 2, 0.2, 0.2, 0.2, 0.2, 3, 3),
    nrow = 4,
    dimnames = list(genes, cts)
  )
  sigma <- array(0, c(4, 4, 2), dimnames = list(genes, genes, cts))
  sigma[,, 1] <- diag(0.05, 4)
  sigma[,, 2] <- diag(0.05, 4)
  sim <- simulate_bulk_from_cells(
    mu,
    sigma,
    p = c(0.8, 0.2),
    n_cells = 10000,
    platform_scale = c(1, 2),
    n_samples = 2,
    seed = 11
  )

  expect_identical(dim(sim$Y), c(4L, 2L))
  expect_equal(unname(colSums(sim$counts)), c(10000, 10000))
  expect_equal(sim$p_realised, sweep(sim$counts, 2, sim$n_cells, "/"))
  expect_equal(sim$mean_scale, sim$platform_scale * sim$n_cells)
  expect_equal(sim$variance_scale, sim$platform_scale^2 * sim$n_cells)
  expect_equal(sim$platform_scale, c(1, 2))
  size_mass <- sim$cell_size * sim$p_realised
  expect_equal(sim$q, sweep(size_mass, 2, colSums(size_mass), "/"))
  expect_true(all(abs(colSums(sim$p_realised) - 1) < 1e-12))
  expect_true(all(abs(colSums(sim$q) - 1) < 1e-12))
})

test_that("the single-cell generator needs more genes than cell types", {
  mu <- matrix(c(1, 0, 0, 1), nrow = 2)
  sigma <- array(rep(diag(2), 2), dim = c(2, 2, 2))
  expect_error(
    simulate_bulk_from_cells(mu, sigma, p = c(0.5, 0.5), n_cells = 10),
    "G > J"
  )
})

test_that("a composition that does not sum to one is rejected", {
  mu <- diag(4)[, 1:2]
  sigma <- array(rep(diag(4), 2), dim = c(4, 4, 2))
  expect_error(
    simulate_bulk_from_cells(mu, sigma, p = c(0.8, 0.8), n_cells = 20),
    "sum to 1"
  )
})
