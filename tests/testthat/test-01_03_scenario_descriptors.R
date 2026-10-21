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
  expect_true("max_cosine" %in% names(out$descriptors))
  expect_true("mean_euclidean" %in% names(out$descriptors))
  expect_true("hellinger_weighted" %in% names(out$descriptors))
  expect_true("riemannian_sigma" %in% names(out$descriptors))
  expect_gt(out$descriptors$mean_euclidean, 0)
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

test_that("describe_simulation_scenario can skip MixSim overlap", {
  theta <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(10, 0), c(0, 10)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  out <- describe_simulation_scenario(theta, include_mixsim = FALSE)
  expect_true(is.na(out$descriptors$mixsim_baromega))
})

test_that("CLD Euclidean gap matches the bivariate small-centroid pair", {
  mu <- cbind(c(20, 22), c(22, 20))
  expect_equal(DeCovarT:::.mean_euclidean_distance(mu), sqrt(8))
  cos_summ <- DeCovarT:::.mean_cosine_summaries(mu)
  expect_gt(cos_summ$max_cosine, 0.99)
})

test_that("Hellinger of two Gaussians is symmetric", {
  mu1 <- c(0, 0)
  mu2 <- c(1, 1)
  s1 <- diag(2)
  s2 <- diag(c(1, 2))
  h12 <- DeCovarT:::.hellinger_gaussian(mu1, mu2, s1, s2)
  h21 <- DeCovarT:::.hellinger_gaussian(mu2, mu1, s2, s1)
  expect_equal(h12, h21)
})


test_that("compute_shannon_entropy: Dirac and uniform", {
  expect_equal(compute_shannon_entropy(c(1, 0, 0)), 0, tolerance = .tol_srr)
  expect_equal(compute_shannon_entropy(rep(1 / 3, 3)), 1, tolerance = .tol_srr)
  expect_equal(compute_shannon_entropy(1), 0, tolerance = .tol_srr)
})

test_that("composition_from_entropy matches a target H-star", {
  p_half <- composition_from_entropy(0.5, n_celltypes = 3L)
  expect_equal(sum(p_half), 1, tolerance = 1e-8)
  expect_equal(compute_shannon_entropy(p_half), 0.5, tolerance = 1e-6)
  expect_equal(
    composition_from_entropy(1, n_celltypes = 3L),
    rep(1 / 3, 3),
    tolerance = 1e-8
  )
  expect_equal(
    composition_from_entropy(0, n_celltypes = 3L),
    c(1, 0, 0)
  )
  p_rare <- composition_from_entropy(0.1, n_celltypes = 3L)
  expect_equal(compute_shannon_entropy(p_rare), 0.1, tolerance = 1e-6)
  expect_gt(p_rare[[1L]], p_rare[[2L]])
})

test_that("compute_shannon_entropy: input checks", {
  expect_error(compute_shannon_entropy(numeric()), "`ratios`")
  expect_error(compute_shannon_entropy(c(0.5, NA)), "`ratios`")
  expect_error(compute_shannon_entropy(c(-0.1, 1.1)), "between 0 and 1")
})
