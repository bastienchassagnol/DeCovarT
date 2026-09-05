test_that("describe_simulation_scenario splits mean and covariance information", {
  mu <- cbind(
    ct1 = c(10, 0),
    ct2 = c(0, 10)
  )
  Sigma <- array(0, dim = c(2, 2, 2))
  Sigma[,, 1] <- diag(c(1, 4))
  Sigma[,, 2] <- diag(c(4, 1))
  theta <- list(
    p = c(0.5, 0.5),
    mu = mu,
    sigma = Sigma
  )
  out <- describe_simulation_scenario(theta)

  expect_named(
    out,
    c("theta_true", "descriptors", "supplementary", "call")
  )
  expect_equal(out$theta_true$p, c(0.5, 0.5))
  expect_type(out$call, "language")
  expect_true("mixsim_baromega" %in% names(out$descriptors))
  expect_true("hellinger" %in% names(out$descriptors))
  expect_true("jeffreys" %in% names(out$supplementary))
  expect_equal(out$descriptors$h_star, 1)
  expect_equal(out$descriptors$n_eff, 2)
  expect_gt(out$descriptors$lambda_min_it, 0)
  expect_gte(out$descriptors$f_cov, 0)
  expect_lte(out$descriptors$f_cov, 1)
  expect_true("kappa_sigma_p" %in% names(out$descriptors))
  expect_false("mixsim_baromega" %in% names(out$supplementary))
  expect_false("hellinger" %in% names(out$supplementary))
})

test_that("identical means remain distinguishable through covariance", {
  Sigma <- array(0, dim = c(2, 2, 2))
  Sigma[,, 1] <- diag(c(1, 4))
  Sigma[,, 2] <- diag(c(4, 1))
  sigma_of_t <- function(t) {
    t^2 * Sigma[,, 1] + (1 - t)^2 * Sigma[,, 2]
  }
  grid <- seq(0, 1, by = 0.05)
  sep <- Inf
  for (a in seq_along(grid)) {
    for (b in seq_along(grid)) {
      if (a != b) {
        sep <- min(
          sep,
          norm(sigma_of_t(grid[[a]]) - sigma_of_t(grid[[b]]), type = "F")
        )
      }
    }
  }
  expect_gt(sep, 0)

  mu <- cbind(ct1 = c(10, 10), ct2 = c(10, 10))
  described <- describe_simulation_scenario(list(
    p = c(0.3, 0.7),
    mu = mu,
    sigma = Sigma
  ))
  expect_gt(described$descriptors$f_cov, 0)
  expect_gt(described$descriptors$f_cov_max, 0.5)
  expect_gt(described$descriptors$lambda_min_it, 0)
})
