# Global performance metrics on the whole composition vector.

test_that("simplex metrics attain the vertex-swap bounds", {
  p <- c(1, 0)
  p_hat <- c(0, 1)
  expect_equal(.tv(p, p_hat), 1)
  expect_equal(.rmse(p, p_hat), 1)
  expect_equal(.maxae(p, p_hat), 1)
  expect_equal(.angular_distance(p, p_hat), 1)
  expect_equal(.sdid(p, p_hat), 1)
  expect_equal(.mae(p, p_hat), 1)
  expect_equal(.tv(p, p), 0)
})

test_that("KKT residual vanishes for a constant ambient score", {
  p <- c(0.2, 0.3, 0.5)
  grad <- c(4, 4, 4)
  expect_equal(.kkt_residual(p, grad), 0, tolerance = 1e-10)
})

test_that("internal metrics match Metrics package on fixed vectors", {
  skip_if_not_installed("Metrics")
  set.seed(1L)
  actual <- c(0.4, 0.6, 0.35, 0.65)
  predicted <- c(0.45, 0.55, 0.4, 0.6)

  expect_equal(.rmse(actual, predicted), Metrics::rmse(actual, predicted))
  expect_equal(.mae(actual, predicted), Metrics::mae(actual, predicted))
})

test_that("internal metrics match Metrics on benchmark-style proportions", {
  skip_if_not_installed("Metrics")
  true_p <- c(0.4, 0.6)
  estimated_p <- c(0.45, 0.55)

  expect_equal(.rmse(true_p, estimated_p), Metrics::rmse(true_p, estimated_p))
  expect_equal(.mae(true_p, estimated_p), Metrics::mae(true_p, estimated_p))
})

test_that("compute_benchmark_metrics uses internal metrics consistently", {
  skip_if_not_installed("Metrics")
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  true_p <- c(0.4, 0.6)
  estimated_p <- c(0.45, 0.55)

  scores <- compute_benchmark_metrics(
    y = drop(mu %*% true_p),
    mean_signature_matrix = mu,
    estimated_p = estimated_p,
    true_ratios = true_p
  )

  expect_equal(
    scores$regression$global$rmse,
    Metrics::rmse(true_p, estimated_p)
  )
  expect_equal(
    scores$regression$global$tv,
    0.5 * sum(abs(true_p - estimated_p))
  )
})
