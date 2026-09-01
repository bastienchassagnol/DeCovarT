# Structured covariance backends must match the dense reference to machine
# precision. The dense path is the ground truth: assemble
# Sigma(p) = sum_j p_j^2 Sigma_j and use base solve / determinant.

dense_reference <- function(p, Sigma) {
  sigma_p <- .compute_global_variance(p, Sigma)
  list(
    matrix = sigma_p,
    logdet = as.numeric(determinant(sigma_p, logarithm = TRUE)$modulus),
    inverse = solve(sigma_p)
  )
}

make_spd <- function(n, seed) {
  withr::with_seed(seed, {
    a <- matrix(stats::rnorm(n * n), n, n)
    crossprod(a) + diag(n)
  })
}

test_that("dense backend matches the base dense reference", {
  Sigma <- array(0, dim = c(3, 3, 2))
  Sigma[,, 1] <- make_spd(3, 1L)
  Sigma[,, 2] <- make_spd(3, 2L)
  p <- c(0.6, 0.4)
  ref <- dense_reference(p, Sigma)

  cov <- new_decovart_covariance(Sigma, "dense")
  r <- c(1, -2, 0.5)

  expect_equal(sigma_logdet(cov, p), ref$logdet)
  expect_equal(sigma_solve(cov, p, r), as.numeric(ref$inverse %*% r))
  expect_equal(
    sigma_quadform(cov, p, r),
    as.numeric(t(r) %*% ref$inverse %*% r)
  )
  expect_equal(
    sigma_trace_precision_times(cov, p, Sigma[,, 1]),
    sum(diag(ref$inverse %*% Sigma[,, 1]))
  )
})

test_that("block backend matches dense for block-diagonal covariance", {
  Sigma <- array(0, dim = c(4, 4, 2))
  Sigma[1:2, 1:2, 1] <- make_spd(2, 11L)
  Sigma[3:4, 3:4, 1] <- make_spd(2, 12L)
  Sigma[1:2, 1:2, 2] <- make_spd(2, 13L)
  Sigma[3:4, 3:4, 2] <- make_spd(2, 14L)
  p <- c(0.7, 0.3)
  ref <- dense_reference(p, Sigma)

  cov <- new_decovart_covariance(Sigma, "block", blocks = c(1, 1, 2, 2))
  r <- c(2, -1, 0.5, 3)

  expect_equal(sigma_logdet(cov, p), ref$logdet)
  expect_equal(sigma_solve(cov, p, r), as.numeric(ref$inverse %*% r))
  expect_equal(
    sigma_quadform(cov, p, r),
    as.numeric(t(r) %*% ref$inverse %*% r)
  )
})

test_that("block backend rejects cross-block covariance", {
  Sigma <- array(0, dim = c(4, 4, 1))
  Sigma[,, 1] <- make_spd(4, 21L)
  expect_error(
    new_decovart_covariance(Sigma, "block", blocks = c(1, 1, 2, 2)),
    "between-block"
  )
})

test_that("band backend matches dense for a tridiagonal covariance", {
  build_band <- function(seed) {
    m <- diag(4) * 3
    off <- withr::with_seed(seed, stats::runif(3, -0.4, 0.4))
    for (i in 1:3) {
      m[i, i + 1] <- off[i]
      m[i + 1, i] <- off[i]
    }
    m
  }
  Sigma <- array(0, dim = c(4, 4, 2))
  Sigma[,, 1] <- build_band(31L)
  Sigma[,, 2] <- build_band(32L)
  p <- c(0.5, 0.5)
  ref <- dense_reference(p, Sigma)

  cov <- new_decovart_covariance(Sigma, "band", bandwidth = 1L)
  r <- c(1, 2, -1, 0.5)

  expect_equal(sigma_logdet(cov, p), ref$logdet, tolerance = 1e-8)
  expect_equal(
    sigma_solve(cov, p, r),
    as.numeric(ref$inverse %*% r),
    tolerance = 1e-8
  )
})

test_that("band backend rejects entries beyond the bandwidth", {
  Sigma <- array(0, dim = c(4, 4, 1))
  Sigma[,, 1] <- make_spd(4, 41L)
  expect_error(
    new_decovart_covariance(Sigma, "band", bandwidth = 1L),
    "bandwidth"
  )
})

test_that("sparse backend matches dense and reuses the symbolic factor", {
  build_sparse <- function(seed) {
    m <- diag(6) * 4
    withr::with_seed(seed, {
      m[1, 4] <- m[4, 1] <- stats::runif(1, -0.3, 0.3)
      m[2, 5] <- m[5, 2] <- stats::runif(1, -0.3, 0.3)
    })
    m
  }
  Sigma <- array(0, dim = c(6, 6, 2))
  Sigma[,, 1] <- build_sparse(51L)
  Sigma[,, 2] <- build_sparse(52L)

  cov <- new_decovart_covariance(Sigma, "sparse")
  r <- c(1, -2, 0.5, 3, -1, 2)

  p1 <- c(0.6, 0.4)
  ref1 <- dense_reference(p1, Sigma)
  expect_equal(sigma_logdet(cov, p1), ref1$logdet, tolerance = 1e-8)
  expect_equal(
    sigma_solve(cov, p1, r),
    as.numeric(ref1$inverse %*% r),
    tolerance = 1e-8
  )

  # Second trial p reuses the cached ordering; result must stay correct.
  p2 <- c(0.3, 0.7)
  ref2 <- dense_reference(p2, Sigma)
  expect_equal(sigma_logdet(cov, p2), ref2$logdet, tolerance = 1e-8)
  expect_equal(
    sigma_solve(cov, p2, r),
    as.numeric(ref2$inverse %*% r),
    tolerance = 1e-8
  )
})

test_that("diag_lowrank backend matches dense via Woodbury", {
  n_genes <- 5L
  rank_r <- 2L
  n_ct <- 2L
  diagonal <- matrix(
    withr::with_seed(61L, stats::runif(n_genes * n_ct, 0.5, 1.5)),
    nrow = n_genes
  )
  loadings <- matrix(
    withr::with_seed(62L, stats::rnorm(n_genes * rank_r)),
    nrow = n_genes
  )
  core <- array(0, dim = c(rank_r, rank_r, n_ct))
  core[,, 1] <- make_spd(rank_r, 63L)
  core[,, 2] <- make_spd(rank_r, 64L)

  cov <- new_decovart_covariance(
    structure = "diag_lowrank",
    diagonal = diagonal,
    loadings = loadings,
    core = core
  )
  Sigma <- cov$Sigma
  p <- c(0.55, 0.45)
  ref <- dense_reference(p, Sigma)
  r <- c(1, -1, 2, 0.5, -0.5)

  expect_equal(sigma_logdet(cov, p), ref$logdet, tolerance = 1e-8)
  expect_equal(
    sigma_solve(cov, p, r),
    as.numeric(ref$inverse %*% r),
    tolerance = 1e-8
  )
  expect_equal(
    sigma_quadform(cov, p, r),
    as.numeric(t(r) %*% ref$inverse %*% r),
    tolerance = 1e-8
  )
})

test_that("diag_lowrank rejects factors that do not reconstruct Sigma", {
  Sigma <- array(0, dim = c(3, 3, 1))
  Sigma[,, 1] <- make_spd(3, 71L)
  expect_error(
    new_decovart_covariance(
      Sigma,
      "diag_lowrank",
      diagonal = matrix(1, 3, 1),
      loadings = matrix(0, 3, 1),
      core = array(1, dim = c(1, 1, 1))
    ),
    "reconstruct"
  )
})

test_that("covariance_structure_from_graph_model maps topologies", {
  expect_identical(
    covariance_structure_from_graph_model("stochastic_block_model"),
    "block"
  )
  expect_identical(
    covariance_structure_from_graph_model("erdos_renyi"),
    "sparse"
  )
  expect_identical(
    covariance_structure_from_graph_model("scale_free"),
    "sparse"
  )
  expect_identical(covariance_structure_from_graph_model("hub"), "diag_lowrank")
  expect_identical(
    covariance_structure_from_graph_model("star"),
    "diag_lowrank"
  )
  expect_identical(covariance_structure_from_graph_model("band"), "band")
})

test_that(".sigma_p_factorisation backend matches the dense factorisation", {
  Sigma <- array(0, dim = c(3, 3, 2))
  Sigma[,, 1] <- make_spd(3, 81L)
  Sigma[,, 2] <- make_spd(3, 82L)
  p <- c(0.4, 0.6)

  dense <- .sigma_p_factorisation(p, Sigma)
  cov <- new_decovart_covariance(Sigma, "dense")
  structured <- .sigma_p_factorisation(p, Sigma, backend = cov)

  expect_equal(structured$log_det, dense$log_det)
  expect_equal(structured$inverse, dense$inverse)
  expect_type(structured$solve, "closure")
  expect_type(structured$quadform, "closure")
})
