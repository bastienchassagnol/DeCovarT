test_that("plot_correlation_Heatmap returns one heatmap per algorithm", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  skip_if_not_installed("viridis")

  metrics <- tibble::tibble(
    correlation_celltype1 = c(0, 0, 0.5, 0.5),
    correlation_celltype2 = c(0, 0.5, 0, 0.5),
    algorithm = "nnls",
    model_mse = c(0.01, 0.02, 0.015, 0.03)
  )

  ht <- plot_correlation_Heatmap(metrics, score_variable = "model_mse")

  expect_type(ht, "list")
  expect_named(ht, "nnls")
  expect_s4_class(ht[["nnls"]], "Heatmap")
})
