sc_inputs <- function(G = 12L, J = 4L) {
  mu <- withr::with_seed(
    1L,
    matrix(
      stats::rexp(G * J, 1 / 50),
      G,
      J,
      dimnames = list(NULL, paste0("ct", seq_len(J)))
    )
  )
  Sigma <- array(0, c(G, G, J))
  for (j in seq_len(J)) {
    Sigma[, , j] <- diag(mu[, j])
  }
  list(mu = mu, Sigma = Sigma)
}

test_that("constructors return likelihood objects, accepted as functions", {
  pop <- population_likelihood()
  sc <- single_cell_likelihood()
  expect_s3_class(pop, "decovart_likelihood")
  expect_identical(c(pop$name, sc$name), c("population", "single_cell"))
  expect_identical(pop$parameters(c("a", "b")), c("a", "b"))
  expect_identical(
    sc$parameters("a"),
    c("a", "mean_scale", "variance_scale", "platform_scale_a", "cell_count_n")
  )
  expect_identical(.as_likelihood(single_cell_likelihood)$name, "single_cell")
  expect_error(.as_likelihood("population"), "likelihood object")
  expect_output(print(sc), "single_cell")
})

test_that("population likelihood scores a fit as the core functions do", {
  mu <- matrix(c(20, 22, 22, 20), 2, dimnames = list(NULL, c("a", "b")))
  Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
  p <- c(a = 0.4, b = 0.6)
  y <- drop(mu %*% p)
  lik <- population_likelihood()
  expect_equal(
    lik$loglik(p, y, mu, Sigma),
    loglik_multivariate(p, y, mu, Sigma)
  )
  expect_equal(
    lik$score(p, y, mu, Sigma),
    gradient_loglik_unconstrained(p, y, mu, Sigma)
  )
})

test_that("every bulk gets the same interval rows, interior or on the boundary", {
  inputs <- sc_inputs()
  mu <- inputs$mu
  Sigma <- inputs$Sigma
  lik <- single_cell_likelihood()
  p_true <- c(ct1 = 0.4, ct2 = 0.3, ct3 = 0.2, ct4 = 0.1)
  # A bulk with noise: a noise-free one has variance scale zero.
  y <- withr::with_seed(2L, {
    V <- 100 * apply(sweep(Sigma, 3, p_true, "*"), 1:2, sum)
    drop(100 * mu %*% p_true + crossprod(chol(V), stats::rnorm(nrow(mu))))
  })
  ci <- list(methods = c("wald", "ilr", "ilr_reml", "oracle"), n_draws = 50L)
  rows <- lik$parameters(names(p_true))
  # Interior, face (one zero), edge (two zeros), vertex (J - 1 zeros).
  estimates <- list(
    interior = p_true,
    face = c(ct1 = 0.45, ct2 = 0.35, ct3 = 0.2, ct4 = 0),
    edge = c(ct1 = 0.6, ct2 = 0.4, ct3 = 0, ct4 = 0),
    vertex = c(ct1 = 1, ct2 = 0, ct3 = 0, ct4 = 0)
  )
  for (k in names(estimates)) {
    out <- .likelihood_intervals(lik, y, estimates[[k]], p_true, mu, Sigma, ci)
    for (method in ci$methods) {
      expect_identical(rownames(out[[method]]), rows, info = paste(k, method))
    }
    finite <- all(is.finite(out$wald))
    expect_identical(finite, k == "interior", info = k)
  }
})

test_that("the coverage table counts boundary estimates", {
  inputs <- sc_inputs()
  mu <- inputs$mu
  Sigma <- inputs$Sigma
  lik <- single_cell_likelihood()
  p_true <- c(ct1 = 0.4, ct2 = 0.3, ct3 = 0.2, ct4 = 0.1)
  # A bulk with noise: a noise-free one has variance scale zero.
  y <- withr::with_seed(2L, {
    V <- 100 * apply(sweep(Sigma, 3, p_true, "*"), 1:2, sum)
    drop(100 * mu %*% p_true + crossprod(chol(V), stats::rnorm(nrow(mu))))
  })
  p_hat <- cbind(p_true, c(0.45, 0.35, 0.2, 0), c(1, 0, 0, 0))
  rownames(p_hat) <- names(p_true)
  ci <- list(methods = "wald")
  intervals <- lapply(seq_len(3L), function(i) {
    .likelihood_intervals(lik, y, p_hat[, i], p_true, mu, Sigma, ci)
  })
  out <- .interval_table(
    intervals,
    p_hat = p_hat,
    p_true = matrix(p_true, 4, 3, dimnames = dimnames(p_hat))
  )
  expect_identical(unique(out$n_boundary), 2L)
  expect_identical(unique(out$n_intervals), 1L)
})
