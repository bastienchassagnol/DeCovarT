##################################################################
##                   bivariate simulations                   ##
##################################################################

test_that("small simulation testing", {
  mu <- matrix(
    c(20, 40, 40, 20),
    nrow = 2,
    dimnames = list(paste0("genes_", 1:2), paste0("cell_type_", 1:2))
  )
  Sigma <- array(
    c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
    dim = c(2, 2, 2),
    dimnames = list(rownames(mu), rownames(mu), colnames(mu))
  )
  simulation_two_genes <- withr::with_seed(
    seed = 3L,
    simulate_bulk_mixture(signature_matrix = mu, Sigma = Sigma, n = 10)
  )
  # Same seed, built by hand: x_j = mu_j + Z chol(Sigma_j), y = (x_1 + x_2) / 2.
  expected <- withr::with_seed(3L, {
    x <- lapply(1:2, function(j) {
      z <- matrix(stats::rnorm(20), 10, 2)
      sweep(z %*% chol(Sigma[, , j]), 2, mu[, j], "+")
    })
    t(0.5 * x[[1]] + 0.5 * x[[2]])
  })
  expect_equal(unname(simulation_two_genes$Y), unname(expected))
})

test_that("simulate_bulk_mixture can truncate negative bulk entries", {
  mu <- matrix(
    c(0.1, 0.1, 0.1, 0.1),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  Sigma <- array(
    c(4, 0, 0, 4, 4, 0, 0, 4),
    dim = c(2, 2, 2),
    dimnames = list(paste0("g", 1:2), paste0("g", 1:2), paste0("ct", 1:2))
  )
  raw <- withr::with_seed(
    1L,
    simulate_bulk_mixture(mu, Sigma, p = c(0.5, 0.5), n = 20L)
  )
  truncated <- withr::with_seed(
    1L,
    simulate_bulk_mixture(
      mu,
      Sigma,
      p = c(0.5, 0.5),
      n = 20L,
      truncate_negative = TRUE
    )
  )
  expect_true(any(raw$Y < 0))
  expect_true(all(truncated$Y >= 0))
})
