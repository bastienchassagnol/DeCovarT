test_that("compute_shannon_entropy: Dirac and uniform", {
  expect_equal(compute_shannon_entropy(c(1, 0, 0)), 0)
  expect_equal(compute_shannon_entropy(rep(1 / 3, 3)), 1)
  expect_equal(compute_shannon_entropy(1), 0)
})

test_that("compute_shannon_entropy: input checks", {
  expect_error(compute_shannon_entropy(numeric()), "`ratios`")
  expect_error(compute_shannon_entropy(c(0.5, NA)), "`ratios`")
  expect_error(compute_shannon_entropy(c(-0.1, 1.1)), "between 0 and 1")
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
  expect_true(is.numeric(ov_far) && length(ov_far) == 1L)
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
  expect_equal(compute_average_jeffreys(theta), 0, tolerance = 1e-10)
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
  expect_equal(j1, j2)
})

test_that("compute_glmnet_gene_scores: returns named G-vector", {
  skip_if_not_installed("glmnet")
  set.seed(42)
  mu <- cbind(
    c(0, 0, 5, 0),
    c(5, 0, 0, 0),
    c(0, 5, 0, 0)
  )
  rownames(mu) <- paste0("g", seq_len(4))
  colnames(mu) <- paste0("ct", seq_len(3))

  scores <- compute_glmnet_gene_scores(mu, n_rep = 15L, nfolds = 5L)
  expect_length(scores, 4L)
  expect_identical(names(scores), rownames(mu))
  expect_true(all(scores >= 0))
  expect_gt(sum(scores), 0)
})

test_that("compute_glmnet_gene_scores: rejects single cell type", {
  skip_if_not_installed("glmnet")
  expect_error(
    compute_glmnet_gene_scores(matrix(1:3, nrow = 3, ncol = 1)),
    "J >= 2"
  )
})
