# Monte Carlo plots of a simulation benchmark.

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
  geoms <- vapply(p_hat$layers, function(ly) class(ly$geom)[[1L]], character(1))
  expect_true(any(grepl("GeomLabel", geoms, fixed = TRUE)))

  p_dots <- plot_mc_raincloud(
    out,
    quantity = "estimate",
    slab = "dotsinterval",
    include_dots = FALSE
  )
  expect_s3_class(p_dots, "ggplot")
  built_dots <- ggplot2::ggplot_build(p_dots)
  expect_type(built_dots$data, "list")
})

test_that("plot_mc_forest shows Wilson coverage whiskers", {
  skip_if_not_installed("nnls")
  skip_if_not_installed("ggplot2")
  out <- .tiny_mc_benchmark()
  p <- plot_mc_forest(
    out,
    facet_rows = "n_genes",
    facet_cols = "cosine",
    metrics = c("bias", "coverage", "mae", "se_sd_ratio")
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
  expect_true(all(p$data$score >= 0 & p$data$score <= 1, na.rm = TRUE))
  expect_true(all(p$data$raw[p$data$metric == "coverage"] >= 0, na.rm = TRUE))
})

test_that(".add_true_ratio_labels puts the abundant type on the right", {
  df <- tibble::tibble(
    cell_type = c("celltype_1", "celltype_2", "celltype_3"),
    p_true = c(0.8, 0.1, 0.1),
    panel = "highly unbalanced"
  )
  p <- .add_true_ratio_labels(
    ggplot2::ggplot(),
    df,
    extra_keys = "panel",
    vlines = TRUE
  )
  lab <- dplyr::bind_rows(lapply(p$layers, function(l) {
    if (inherits(l$geom, "GeomLabel")) l$data else NULL
  }))
  expect_identical(nrow(lab), 3L)
  ct1 <- lab[lab$cell_type == "celltype_1", ]
  expect_true(ct1$lab_x > ct1$p_true)
  expect_identical(ct1$hjust, 0)
  others <- lab[lab$cell_type != "celltype_1", ]
  expect_true(all(others$lab_x < others$p_true))
  expect_true(all(others$hjust == 1))
  expect_identical(
    sum(vapply(
      p$layers,
      function(l) inherits(l$geom, "GeomVline"),
      logical(1)
    )),
    3L
  )
})
