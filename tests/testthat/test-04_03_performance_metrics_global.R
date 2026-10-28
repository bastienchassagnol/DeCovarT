# Global performance metrics on the whole composition vector.

.global_scores <- function(true_p, estimated_p, y = NULL) {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  if (is.null(y)) {
    y <- drop(mu %*% c(0.4, 0.6))
  }
  compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = estimated_p,
    true_ratios = true_p
  )$regression$global
}

test_that("simplex metrics attain the vertex-swap bounds", {
  swap <- .global_scores(c(1, 0), c(0, 1))
  expect_equal(swap$tv, 1)
  expect_equal(swap$rmse, 1)
  expect_equal(swap$maxae, 1)
  expect_equal(swap$angular, 1)
  expect_equal(swap$sdid, 1)
  exact <- .global_scores(c(1, 0), c(1, 0))
  expect_equal(exact$tv, 0)
  expect_equal(exact$angular, 0)
  expect_equal(exact$sdid, 0)
})

test_that("global scores match hand computations", {
  true_p <- c(0.4, 0.6)
  estimated_p <- c(0.45, 0.55)
  scores <- .global_scores(true_p, estimated_p)
  cos_angle <- sum(true_p * estimated_p) /
    sqrt(sum(true_p^2) * sum(estimated_p^2))
  expect_equal(scores$tv, 0.05)
  expect_equal(scores$rmse, 0.05)
  expect_equal(scores$maxae, 0.05)
  expect_equal(scores$angular, (2 / pi) * acos(cos_angle))
  expect_equal(scores$sdid, sqrt(sqrt(1 - cos_angle^2)))
})

test_that("reconstitution scores compare mu p_hat with the bulk", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  y <- c(21, 21.5)
  scores <- .global_scores(NULL, c(0.45, 0.55), y = y)
  error <- y - drop(mu %*% c(0.45, 0.55))
  expect_equal(scores$rmse, sqrt(mean(error^2)))
  expect_equal(scores$reconstitution_mae, mean(abs(error)))
  expect_true(is.na(scores$tv))
})

test_that("internal RMSE and MAE match the Metrics package", {
  skip_if_not_installed("Metrics")
  true_p <- c(0.4, 0.6)
  estimated_p <- c(0.45, 0.55)
  scores <- .global_scores(true_p, estimated_p)
  expect_equal(scores$rmse, Metrics::rmse(true_p, estimated_p))
  expect_equal(scores$tv, Metrics::mae(true_p, estimated_p))
})

test_that("KKT residual vanishes for a constant ambient score", {
  p <- c(0.2, 0.3, 0.5)
  grad <- c(4, 4, 4)
  expect_equal(kkt_residual(p, grad), 0, tolerance = 1e-10)
})
