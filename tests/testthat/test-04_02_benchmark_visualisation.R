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

.tiny_mc_benchmark <- function() {
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  theta <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
  scenario_config <- tibble::tibble(
    n_genes = c(2L, 2L),
    cosine = c(0, 0.5),
    true_theta = list(theta, theta)
  )
  withr::with_seed(
    11L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls)
      ),
      n = 6L,
      cores = 1L
    )
  )
}

test_that("pivot_mc_estimates stacks Monte Carlo proportion draws", {
  skip_if_not_installed("nnls")
  out <- .tiny_mc_benchmark()
  long <- pivot_mc_estimates(out)
  expect_s3_class(long, "tbl_df")
  needed <- c(
    "algorithm",
    "cell_type",
    "estimate",
    "p_true",
    "error",
    "n_genes",
    "cosine"
  )
  expect_true(all(needed %in% names(long)))
  n_scen <- nrow(out$config)
  n_opt <- nrow(out$optimisation)
  n_ct <- length(unique(long$cell_type))
  expect_equal(n_opt, n_scen * 6L)
  expect_equal(nrow(long), n_opt * n_ct)
  expect_equal(long$error, long$estimate - long$p_true)
  expect_true(all(long$algorithm == "nnls"))
})

test_that("plot_mc_raincloud builds a horizontal raincloud ggplot", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdist")
  out <- .tiny_mc_benchmark()
  p_err <- plot_mc_raincloud(
    out,
    quantity = "error",
    facet_rows = "n_genes",
    facet_cols = "cosine"
  )
  expect_s3_class(p_err, "ggplot")
  built_err <- ggplot2::ggplot_build(p_err)
  expect_type(built_err$data, "list")
  expect_gt(length(built_err$data), 0L)
  expect_true("FacetGrid" %in% class(p_err$facet))

  p_hat <- plot_mc_raincloud(
    out,
    quantity = "estimate",
    facet_cols = "cosine"
  )
  expect_s3_class(p_hat, "ggplot")
  built_hat <- ggplot2::ggplot_build(p_hat)
  expect_type(built_hat$data, "list")
})

test_that("plot_mc_forest shows Wilson coverage whiskers", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("ggplot2")
  out <- .tiny_mc_benchmark()
  p <- plot_mc_forest(
    out,
    facet_rows = "n_genes",
    facet_cols = "cosine",
    metrics = c("bias", "coverage", "mae")
  )
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_type(built$data, "list")
  expect_true("FacetGrid" %in% class(p$facet))
  mc <- out$monte_carlo
  expect_true(all(c("coverage_lower", "coverage_upper") %in% names(mc)))
  expect_true(all(mc$coverage_interval == "wilson"))
  expect_true(all(mc$coverage_lower <= mc$coverage, na.rm = TRUE))
  expect_true(all(mc$coverage <= mc$coverage_upper, na.rm = TRUE))
})

test_that(".check_plot_dependencies succeeds when ggplot2 is present", {
  skip_if_not_installed("ggplot2")
  expect_true(DeCovarT:::.check_plot_dependencies(need_ggdist = FALSE))
})

test_that(".check_plot_dependencies succeeds when ggdist is present", {
  skip_if_not_installed("ggdist")
  expect_true(DeCovarT:::.check_plot_dependencies(need_ggdist = TRUE))
})

.tiny_mc_benchmark_two_algos <- function() {
  genes <- paste0("g", 1:2)
  cts <- paste0("ct", 1:2)
  mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
  Sigma <- array(
    c(1, 0, 0, 1, 1, 0, 0, 1),
    dim = c(2, 2, 2),
    dimnames = list(genes, genes, cts)
  )
  theta <- list(p = c(0.5, 0.5), mu = mu, sigma = Sigma)
  scenario_config <- tibble::tibble(
    n_genes = 2L,
    cosine = 0,
    true_theta = list(theta)
  )
  withr::with_seed(
    17L,
    run_simulation_benchmark(
      scenario_config = scenario_config,
      deconvolution_functions = list(
        "nnls" = list(FUN = deconvolute_ratios_nnls),
        "rlm" = list(FUN = deconvolute_ratios_rlm)
      ),
      n = 6L,
      cores = 1L
    )
  )
}

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
  out <- .tiny_mc_benchmark_two_algos()
  p <- plot_algorithm_similarity(out, dendrogram = TRUE)
  dend <- attr(p, "dendrogram")
  expect_s3_class(dend, "ggplot")
})

test_that("plot_mc_metric_dots uses one colour scale in [0, 1]", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("ggplot2")
  out <- .tiny_mc_benchmark_two_algos()
  p <- plot_mc_metric_dots(
    out,
    metrics = c("rmse", "coverage")
  )
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_type(built$data, "list")
  expect_true("FacetGrid" %in% class(p$facet))
  scored <- DeCovarT:::.mc_metric_dot_table(out, c("rmse", "coverage"))
  expect_true(all(scored$raw[scored$metric == "coverage"] >= 0, na.rm = TRUE))
})
