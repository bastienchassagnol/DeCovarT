# Figure books of scripts/fig02_bivariate_toy.R,
# scripts/fig03_variance_covariance_driven.R and
# scripts/fig05_single_cell_hierarchical.R.
#
#   plot_data <- figure_plot_data(read_benchmark(file), "fig03")
#   render_figure_books(plot_data, "output/fig03")
#   render_figure_from_data("ci_width_forest", plot_data, "output/fig03")
#
# `plot_data` is a list keyed by file stem; each entry holds the plot
# function, its data and the page size, so any book is redrawn from
# figXX_plot_data.rds alone (after sourcing this file).

source(file.path("scripts", "benchmark_metrics", "outcome_metrics.R"))
source(file.path("scripts", "visualisation", "visualisation_helpers.R"))
source(file.path("scripts", "visualisation", "density_books.R"))
source(file.path("scripts", "visualisation", "performance_books.R"))
source(file.path("scripts", "visualisation", "scenario_books.R"))

.book <- function(builder, data, folder, width, height, ...) {
  list(
    builder = builder,
    data = data,
    arguments = list(...),
    folder = folder,
    width = width,
    height = height
  )
}

figure_plot_data <- function(benchmark, figure = c("fig02", "fig03", "fig05")) {
  figure <- match.arg(figure)
  layout <- .scenario_layout(benchmark$config, figure)
  density <- "density_visualisations"
  performance <- "performance_visualisations"
  books <- list()
  if (identical(figure, "fig02")) {
    densities <- .data_gene_densities(benchmark, layout)
    profiles <- .data_loglik_profiles(benchmark, layout)
    books$purified_density <- .book(
      "plot_gene_densities",
      densities,
      density,
      11,
      11,
      which = "purified"
    )
    books$bulk_density <- .book(
      "plot_gene_densities",
      densities,
      density,
      11,
      11,
      which = "bulk"
    )
    books$loglik_ilr_profile <- .book(
      "plot_loglik_profiles",
      profiles,
      density,
      11,
      9,
      surface = "observed"
    )
    books$loglik_expected_ilr_profile <- .book(
      "plot_loglik_profiles",
      profiles,
      density,
      11,
      9,
      surface = "expected"
    )
  } else {
    books$mean_signature <- .book(
      "plot_mean_signature",
      .data_mean_signature(benchmark, layout, figure),
      density,
      7,
      10
    )
    books$network_topologies <- .book(
      "plot_network_topologies",
      .data_network_topologies(benchmark, layout, figure),
      density,
      13,
      11
    )
    metrics <- .data_scenario_metrics(benchmark, layout, figure)
    books$scenario_metrics <- .book(
      "plot_scenario_metrics",
      metrics,
      density,
      NA,
      NA
    )
    books$scenario_metrics_table <- .book(
      "write_scenario_metrics_table",
      metrics,
      density,
      NA,
      NA
    )
    books$latent_projections <- .book(
      "plot_latent_projections",
      .data_latent_projections(benchmark, layout, figure),
      density,
      12,
      6.5
    )
    books$pairwise_network_aware_distances <- .book(
      "plot_pairwise_distances",
      .data_pairwise_distances(benchmark, layout, figure),
      density,
      10,
      if (identical(figure, "fig05")) 14 else 9,
      figure = figure
    )
  }
  surfaces <- .data_loglik_surfaces(
    benchmark,
    layout,
    figure,
    n_grid = if (identical(figure, "fig05")) 35L else 60L
  )
  books$loglik_surface_p <- .book(
    "plot_loglik_surfaces",
    surfaces,
    density,
    11,
    11.5,
    surface = "observed"
  )
  books$loglik_expected_surface_p <- .book(
    "plot_loglik_surfaces",
    surfaces,
    density,
    11,
    11.5,
    surface = "expected"
  )
  books$loglik_rgl <- .book(
    "plot_loglik_rgl",
    surfaces[c("surface", "markers", "axis_labels")],
    density,
    11,
    10
  )

  errors <- .data_errors(benchmark, layout)
  whitened <- .data_whitened(benchmark, layout, figure)
  estimates <- .data_estimates(benchmark, layout, figure)
  optimisation <- .data_optimisation(benchmark, layout)
  # CIBERSORT's error in the single-cell scenario is an order of magnitude
  # above the others and would squeeze the shared accuracy scales.
  exclude <- if (identical(figure, "fig05")) "CIBERSORT" else character()
  books$heatmap_rmse <- .book(
    "plot_metric_tiles",
    errors,
    performance,
    14,
    10,
    metric = "rmse",
    exclude = exclude
  )
  books$heatmap_aitchison <- .book(
    "plot_metric_tiles",
    errors,
    performance,
    14,
    10,
    metric = "aitchison",
    exclude = exclude
  )
  books$rmse_colour_aitchison_size_dots <- .book(
    "plot_rmse_aitchison_dots",
    errors,
    performance,
    14,
    10,
    exclude = exclude
  )
  books$ci_width_forest <- .book(
    "plot_ci_forest",
    .data_ci_forest(benchmark, layout, figure),
    performance,
    12,
    8
  )
  books$p_hat_monte_carlo_raincloud <- .book(
    "plot_estimate_laws",
    estimates,
    performance,
    12,
    9,
    style = "raincloud"
  )
  books$p_hat_monte_carlo_quantile_dots <- .book(
    "plot_estimate_laws",
    estimates,
    performance,
    12,
    9,
    style = "quantile_dots"
  )
  books$wald_normality_qq <- .book(
    "plot_wald_normality_qq",
    whitened,
    performance,
    12,
    9
  )
  books$wald_mahalanobis_chi2_qq <- .book(
    "plot_wald_chi2_qq",
    whitened,
    performance,
    12,
    9
  )
  books$estimate_agreement_hellinger_ward <- .book(
    "plot_estimate_agreement",
    .data_agreement(benchmark, layout),
    performance,
    10,
    11
  )
  books$kkt_residual_raincloud <- .book(
    "plot_kkt_raincloud",
    optimisation,
    performance,
    12,
    9
  )
  books$convergence_theoretical_vs_numerical_failures <- .book(
    "plot_convergence",
    optimisation,
    performance,
    12,
    9
  )
  books$runtime <- .book(
    "plot_resource",
    optimisation,
    performance,
    12,
    9,
    metric = "elapsed_sec"
  )
  books$memory <- .book(
    "plot_resource",
    optimisation,
    performance,
    12,
    9,
    metric = "memory_mib"
  )
  if (identical(figure, "fig05")) {
    books$ci_width_forest_nuisance <- .book(
      "plot_ci_forest_nuisance",
      .data_ci_forest_nuisance(benchmark, layout),
      performance,
      13,
      10
    )
    books$nuisance_scale_ratio <- .book(
      "plot_nuisance_ratio",
      .data_nuisance(benchmark, layout),
      performance,
      11,
      10
    )
  }
  books
}

render_figure_from_data <- function(stem, plot_data, out_dir) {
  book <- plot_data[[stem]]
  if (is.null(book)) {
    stop("No book `", stem, "` in the plot data.", call. = FALSE)
  }
  if (startsWith(book$builder, "write_")) {
    # Tables write their own files (HTML and PDF) next to the books.
    return(do.call(
      book$builder,
      c(
        list(book$data, file.path(out_dir, book$folder, stem)),
        book$arguments
      )
    ))
  }
  pages <- do.call(book$builder, c(list(book$data), book$arguments))
  width <- book$width
  height <- book$height
  if (is.na(width)) {
    width <- pages[[1L]]$width
    height <- pages[[1L]]$height
  }
  .save_pages(
    pages,
    file.path(out_dir, book$folder, paste0(stem, ".pdf")),
    width,
    height
  )
}

render_figure_books <- function(plot_data, out_dir) {
  for (stem in names(plot_data)) {
    render_figure_from_data(stem, plot_data, out_dir)
  }
  invisible(out_dir)
}
