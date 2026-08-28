test_that("internal metrics match Metrics package on fixed vectors", {
  skip_if_not_installed("Metrics")
  set.seed(1L)
  actual <- c(0.4, 0.6, 0.35, 0.65)
  predicted <- c(0.45, 0.55, 0.4, 0.6)

  expect_equal(.mse(actual, predicted), Metrics::mse(actual, predicted))
  expect_equal(.rmse(actual, predicted), Metrics::rmse(actual, predicted))
  expect_equal(.mae(actual, predicted), Metrics::mae(actual, predicted))
  expect_equal(.rse(actual, predicted), Metrics::rse(actual, predicted))
})

test_that("internal metrics match Metrics on benchmark-style proportions", {
  skip_if_not_installed("Metrics")
  true_p <- c(0.4, 0.6)
  estimated_p <- c(0.45, 0.55)

  expect_equal(.mse(true_p, estimated_p), Metrics::mse(true_p, estimated_p))
  expect_equal(.rmse(true_p, estimated_p), Metrics::rmse(true_p, estimated_p))
  expect_equal(.mae(true_p, estimated_p), Metrics::mae(true_p, estimated_p))
  expect_equal(.rse(true_p, estimated_p), Metrics::rse(true_p, estimated_p))
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

  expect_equal(scores$model_mse, Metrics::mse(true_p, estimated_p))
  expect_equal(scores$model_rmse, Metrics::rmse(true_p, estimated_p))
  expect_equal(scores$model_mae, Metrics::mae(true_p, estimated_p))
  expect_equal(
    scores$model_coef_determination,
    max(0, 1 - Metrics::rse(true_p, estimated_p))
  )
})
