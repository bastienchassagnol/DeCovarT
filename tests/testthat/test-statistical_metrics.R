test_that("compute_shannon_entropy: Dirac and uniform", {
  expect_equal(compute_shannon_entropy(c(1, 0, 0)), 0, tolerance = .tol_srr)
  expect_equal(compute_shannon_entropy(rep(1 / 3, 3)), 1, tolerance = .tol_srr)
  expect_equal(compute_shannon_entropy(1), 0, tolerance = .tol_srr)
})

test_that("compute_shannon_entropy: input checks", {
  expect_error(compute_shannon_entropy(numeric()), "`ratios`")
  expect_error(compute_shannon_entropy(c(0.5, NA)), "`ratios`")
  expect_error(compute_shannon_entropy(c(-0.1, 1.1)), "between 0 and 1")
})

test_that("check_true_theta: accepts valid theta", {
  theta <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(3, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_true(check_true_theta(theta))
  expect_true(check_true_theta(theta, second_moment = "SIGMA"))
})

test_that("check_true_theta: accepts J x N proportions", {
  theta <- list(
    p = cbind(c(0.6, 0.4), c(0.4, 0.6)),
    mu = cbind(c(0, 0), c(3, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_true(check_true_theta(theta))
  ov <- compute_average_overlap(theta)
  expect_type(ov, "double")
  expect_length(ov, 1L)
})

test_that("repair_simplex renormalises and rejects invalid inputs", {
  p_near <- c(0.2, 0.3, 0.5 + 1e-12)
  expect_equal(
    repair_simplex(p_near),
    p_near / sum(p_near),
    tolerance = .tol_srr
  )
  expect_equal(sum(repair_simplex(c(2, 2, 0))), 1, tolerance = .tol_srr)
  expect_error(repair_simplex(numeric()), "`p`")
  expect_error(repair_simplex(c(0.5, NA)), "`p`")
  expect_error(repair_simplex(c(-1, 2)), "negative")
  expect_error(repair_simplex(c(0, 0)), "positive")
})

test_that("check_true_theta: rejects bad dims", {
  expect_error(
    check_true_theta(list(mu = matrix(1, 2, 2), p = c(0.5, 0.5))),
    "sigma"
  )
})

test_that("compute_average_overlap: closer means overlap more", {
  set.seed(1)
  far <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(8, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  near <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(1, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )

  ov_far <- compute_average_overlap(far)
  ov_near <- compute_average_overlap(near)
  expect_type(ov_far, "double")
  expect_length(ov_far, 1L)
  expect_gte(ov_far, 0)
  expect_lte(ov_far, 1)
  expect_lt(ov_far, ov_near)
})

test_that("compute_average_overlap: J mismatch errors", {
  theta <- list(
    p = c(0.5, 0.5),
    mu = cbind(c(0, 0), c(1, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_error(compute_average_overlap(theta, J = 3L), "`J`")
})

test_that("compute_average_jeffreys: identical means give ~0", {
  theta <- list(
    mu = cbind(c(0, 0), c(0, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_equal(compute_average_jeffreys(theta), 0, tolerance = .tol_srr)
})

test_that("compute_average_jeffreys: mean shift increases divergence", {
  close <- list(
    mu = cbind(c(0, 0), c(0.5, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  far <- list(
    mu = cbind(c(0, 0), c(4, 0)),
    sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  )
  expect_lt(compute_average_jeffreys(close), compute_average_jeffreys(far))
})

test_that("compute_average_jeffreys: equi-balanced default when p omitted", {
  theta <- list(
    mu = cbind(c(0, 0), c(2, 0), c(0, 2)),
    sigma = array(c(diag(2), diag(2), diag(2)), dim = c(2, 2, 3))
  )
  j1 <- compute_average_jeffreys(theta)
  theta$p <- rep(1 / 3, 3)
  j2 <- compute_average_jeffreys(theta)
  expect_equal(j1, j2, tolerance = .tol_srr)
})

test_that("compute_glmnet_gene_scores: returns named G-vector", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  G <- 4L
  J <- 3L
  N <- 15L
  profiles <- array(0, dim = c(G, J, N))
  dimnames(profiles) <- list(
    paste0("g", seq_len(G)),
    paste0("ct", seq_len(J)),
    NULL
  )
  for (j in seq_len(J)) {
    profiles[j, j, ] <- 5 + stats::rnorm(N, sd = 0.25)
  }
  labels <- paste0("ct", seq_len(J))

  scores <- compute_glmnet_gene_scores(profiles, labels)
  expect_length(scores, G)
  expect_named(scores, dimnames(profiles)[[1L]])
  expect_true(all(scores >= 0))
  expect_gt(sum(scores), 0)
})

test_that("compute_glmnet_gene_scores: rejects single cell type", {
  skip_if_not_installed("glmnet")
  profiles <- array(1, dim = c(3L, 1L, 5L))
  expect_error(
    compute_glmnet_gene_scores(profiles, "ct1"),
    "J >= 2"
  )
})
