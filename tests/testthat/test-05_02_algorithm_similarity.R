# Algorithm-similarity correlations and heatmap.

test_that("algorithm_similarity is a square Pearson matrix in long form", {
  skip_if_not_installed("nnls")
  out <- .tiny_mc_benchmark_two_algos()
  sim <- algorithm_similarity(out)
  expect_s3_class(sim, "tbl_df")
  expect_true(all(
    c("algorithm_x", "algorithm_y", "correlation") %in% names(sim)
  ))
  expect_equal(nrow(sim), 4L)
  diag_r <- sim$correlation[sim$algorithm_x == sim$algorithm_y]
  expect_equal(diag_r, rep(1, length(diag_r)))
  expect_error(algorithm_similarity(.tiny_mc_benchmark()), "two algorithms")
})

test_that("plot_algorithm_similarity is a clustered geom_tile ggplot", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("ggplot2")
  out <- .tiny_mc_benchmark_two_algos()
  p <- plot_algorithm_similarity(out)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_type(built$data, "list")
})

test_that("plot_algorithm_similarity can attach a ggdendro dendrogram", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdendro")
  skip_if_not_installed("cowplot")
  out <- .tiny_mc_benchmark_two_algos()
  p <- plot_algorithm_similarity(out, dendrogram = TRUE)
  dend <- attr(p, "dendrogram")
  expect_s3_class(dend, "ggplot")
})
