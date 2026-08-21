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

test_that("plot_correlation_Heatmap writes a PDF via tempfile (G4.0)", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  skip_if_not_installed("viridis")

  metrics <- tibble::tibble(
    correlation_celltype1 = c(0, 0, 0.5, 0.5),
    correlation_celltype2 = c(0, 0.5, 0, 0.5),
    algorithm = "nnls",
    model_mse = c(0.01, 0.02, 0.015, 0.03)
  )

  withr::with_tempfile("tf", fileext = ".pdf", {
    ht <- plot_correlation_Heatmap(
      metrics,
      score_variable = "model_mse",
      file = tf
    )
    expect_true(file.exists(tf))
    expect_gt(file.info(tf)$size, 0)
    expect_s4_class(ht[["nnls"]], "Heatmap")
  })
})

test_that(".check_heatmap_dependencies succeeds when Suggests present", {
  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("circlize")
  skip_if_not_installed("viridis")
  expect_true(DeCovarT:::.check_heatmap_dependencies())
})
