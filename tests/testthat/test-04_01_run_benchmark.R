test_that("compute_benchmark_metrics scores against true ratios", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  true_p <- c(0.4, 0.6)
  y <- drop(mu %*% true_p)
  estimated_p <- c(0.45, 0.55)

  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = estimated_p,
    true_ratios = true_p
  )

  expect_named(
    scores,
    c(
      "model_mse",
      "model_rmse",
      "model_mae",
      "model_coef_determination",
      "model_coef_determination_adjusted",
      "model_cor"
    )
  )
  expect_gt(scores$model_mse, 0)
  expect_lte(scores$model_mse, Metrics::mse(true_p, true_p) + 1)
})

test_that("compute_benchmark_metrics scores reconstituted bulk without truth", {
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(paste0("g", 1:2), paste0("ct", 1:2))
  )
  estimated_p <- c(0.5, 0.5)
  y <- drop(mu %*% c(0.4, 0.6))

  scores <- compute_benchmark_metrics(
    y = y,
    mean_signature_matrix = mu,
    estimated_p = estimated_p,
    true_ratios = NULL
  )

  expect_named(
    scores,
    c("model_mse", "model_rmse", "model_mae", "model_cor")
  )
  expect_true(is.finite(scores$model_mse))
})

test_that("deconvolute_ratios returns a tidy table for nnls", {
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(
    c(20, 22, 22, 20),
    nrow = 2,
    dimnames = list(genes, cts)
  )
  bulk <- cbind(
    drop(mu %*% c(0.5, 0.5)),
    drop(mu %*% c(0.3, 0.7))
  )
  dimnames(bulk) <- list(genes, paste0("sample_", 1:2))

  out <- deconvolute_ratios(
    signature_matrix = mu,
    bulk_expression = bulk,
    deconvolution_functions = list(
      "nnls" = list(FUN = deconvolute_ratios_nnls)
    ),
    cores = 1
  )

  expect_true(tibble::is_tibble(out))
  expect_identical(nrow(out), 2L)
  expect_true("algorithm" %in% names(out))
  expect_true(all(c("ct1", "ct2") %in% names(out)))
})
