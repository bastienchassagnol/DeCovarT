# Figure books of scripts/fig02_bivariate_toy.R,
# scripts/fig03_variance_covariance_driven.R,
# scripts/fig05_bivariate_single_cell.R and
# scripts/fig06_single_cell_hierarchical.R.
#
#   plot_data <- figure_plot_data(readRDS(file), "fig03")
#   render_figure_books(plot_data, "output/fig03")
#   render_figure_from_data("ci_width_forest", plot_data, "output/fig03")
#
# `plot_data` is a list keyed by file stem; each entry holds the plot
# function, its data and the page size, so any book is redrawn from
# figXX_plot_data.rds alone (after sourcing this file). fig05 is split into
# one entry per page, keyed `<stem>/<page>`, drawn as one SVG with its RDS
# (the entry) next to it.

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

figure_plot_data <- function(
  benchmark,
  figure = c("fig02", "fig03", "fig05", "fig06")
) {
  figure <- match.arg(figure)
  layout <- .scenario_layout(benchmark$config, figure)
  density <- "density_visualisations"
  performance <- "performance_visualisations"
  books <- list()
  if (figure %in% c("fig02", "fig05")) {
    densities <- .data_gene_densities(benchmark, layout, figure = figure)
    profiles <- .data_loglik_profiles(benchmark, layout, figure = figure)
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
      if (identical(figure, "fig06")) 14 else 9,
      figure = figure
    )
  }
  surfaces <- .data_loglik_surfaces(
    benchmark,
    layout,
    figure,
    n_grid = if (identical(figure, "fig06")) 35L else 60L
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
  if (!identical(figure, "fig05")) {
    books$loglik_rgl <- .book(
      "plot_loglik_rgl",
      surfaces[c("surface", "markers", "axis_labels")],
      density,
      11,
      10
    )
  }

  errors <- .data_errors(benchmark, layout)
  # fig05 pages are drawn one at a time: their colour scales are fixed here.
  shared <- function(x) if (identical(figure, "fig05")) x
  whitened <- .data_whitened(benchmark, layout, figure)
  estimates <- .data_estimates(benchmark, layout, figure)
  optimisation <- .data_optimisation(benchmark, layout)
  # These methods' errors in the single-cell scenario are far above the
  # others and would squeeze the shared accuracy scales. They are removed
  # from the data, including their factor levels, so no panel remains.
  exclude <- if (identical(figure, "fig06")) {
    c("CIBERSORT", "lsei", "MEAD")
  } else {
    character()
  }
  books$heatmap_rmse <- .book(
    "plot_metric_tiles",
    errors,
    performance,
    14,
    10,
    metric = "rmse",
    exclude = exclude,
    limits = shared(.positive_range(errors$rmse))
  )
  books$heatmap_aitchison <- .book(
    "plot_metric_tiles",
    errors,
    performance,
    14,
    10,
    metric = "aitchison",
    exclude = exclude,
    limits = shared(.positive_range(errors$aitchison))
  )
  books$rmse_colour_aitchison_size_dots <- .book(
    "plot_rmse_aitchison_dots",
    errors,
    performance,
    14,
    10,
    exclude = exclude,
    rmse_limits = shared(.positive_range(errors$rmse)),
    size_limits = shared(range(errors$aitchison, finite = TRUE))
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
    books$loglik_mean_shift_slice <- .book(
      "plot_loglik_slices",
      .data_loglik_slices(benchmark, "mean_shift"),
      density,
      11,
      12
    )
    books$loglik_correlation_slice <- .book(
      "plot_loglik_slices",
      .data_loglik_slices(benchmark, "correlation"),
      density,
      11,
      10
    )
    books$boundary_lrt <- .book(
      "plot_boundary_lrt",
      .data_boundary_lrt(benchmark, layout),
      performance,
      11,
      11
    )
    mechanism <- .data_rmse_mechanism(benchmark, layout, errors)
    books$heatmap_rmse_ratio_lsei <- .book(
      "plot_rmse_ratio_tiles",
      mechanism[c("errors", "paired")],
      performance,
      17,
      6.5,
      which = "lsei"
    )
    books$heatmap_rmse_ratio_population <- .book(
      "plot_rmse_ratio_tiles",
      mechanism[c("paired", "errors")],
      performance,
      10,
      7,
      which = "population"
    )
    books$rmse_overlap_scatter <- .book(
      "plot_rmse_overlap_scatter",
      mechanism["scatter"],
      performance,
      16,
      9
    )
    # Drawn page by page: one SVG per page (local) and one PDF per stem.
    books <- lapply(books, function(book) {
      book$format <- "pages"
      book
    })
  }
  if (identical(figure, "fig06")) {
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

# Pages of one book: every table that has a `page` column is cut to one
# page (and, for the interval forests, to one scenario). Returns a list of
# books named by page slug (the book's own stem when it has no pages), with
# a `pages` attribute (`page_id`, `page`, `ID`) to find a page by name.
.split_book_by_page <- function(book, stem) {
  slug <- function(x) {
    x <- gsub(",", "", tolower(as.character(x)), fixed = TRUE)
    x <- sub(" cells per bulk$", "", x)
    gsub("^_|_$", "", gsub("[^a-z0-9]+", "_", x))
  }
  tables <- if (is.data.frame(book$data)) list(book$data) else book$data
  has_page <- vapply(
    tables,
    function(t) is.data.frame(t) && "page" %in% names(t),
    logical(1)
  )
  if (!any(has_page)) {
    out <- stats::setNames(list(book), stem)
    attr(out, "pages") <- data.frame(page_id = stem, page = NA, ID = NA)
    return(out)
  }
  by <- if (identical(stem, "ci_width_forest")) c("page", "ID") else "page"
  keys <- unique(as.data.frame(tables[[which(has_page)[[1L]]]])[by])
  out <- list()
  for (k in seq_len(nrow(keys))) {
    cut <- function(t) {
      if (!is.data.frame(t) || !all(by %in% names(t))) {
        return(t)
      }
      keep <- Reduce(`&`, lapply(by, function(b) t[[b]] == keys[[b]][[k]]))
      t[keep, , drop = FALSE]
    }
    page_book <- book
    page_book$data <- if (is.data.frame(book$data)) {
      cut(book$data)
    } else {
      lapply(book$data, cut)
    }
    key <- paste(
      vapply(by, function(b) slug(keys[[b]][[k]]), ""),
      collapse = "_"
    )
    out[[key]] <- page_book
  }
  attr(out, "pages") <- data.frame(
    page_id = names(out),
    page = as.character(keys$page),
    ID = if ("ID" %in% by) keys$ID else NA
  )
  out
}

# Draw one book. A "pages" book (fig05) writes every page to
# `<folder>/<stem>/<page_id>.svg` (kept locally), all pages to one
# `<folder>/<stem>.pdf`, and the book itself, in long format with its page
# table, to `<folder>/<stem>.rds`; `page` redraws a single page (SVG only).
render_figure_from_data <- function(stem, plot_data, out_dir, page = NULL) {
  book <- if (is.list(plot_data) && is.null(plot_data$builder)) {
    plot_data[[stem]]
  } else {
    plot_data
  }
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
  draw <- function(b) {
    pages <- Filter(
      Negate(is.null),
      do.call(b$builder, c(list(b$data), b$arguments))
    )
    width <- b$width
    height <- b$height
    if (is.na(width)) {
      width <- pages[[1L]]$width
      height <- pages[[1L]]$height
    }
    list(pages = pages, width = width, height = height)
  }
  if (identical(book$format, "pages")) {
    parts <- .split_book_by_page(book, stem)
    index <- attr(parts, "pages")
    if (!is.null(page)) {
      if (!page %in% names(parts)) {
        stop("No page `", page, "` in `", stem, "`.", call. = FALSE)
      }
      parts <- parts[page]
    }
    base <- file.path(out_dir, book$folder, stem)
    all_pages <- list()
    files <- character()
    for (id in names(parts)) {
      d <- draw(parts[[id]])
      f <- if (length(d$pages) == 1L) {
        file.path(base, paste0(id, ".svg"))
      } else {
        file.path(base, sprintf("%s_%d.svg", id, seq_along(d$pages)))
      }
      for (k in seq_along(d$pages)) {
        .save_pages(d$pages[k], f[[k]], d$width, d$height)
      }
      all_pages <- c(all_pages, d$pages)
      files <- c(files, f)
    }
    if (is.null(page)) {
      .save_pages(all_pages, paste0(base, ".pdf"), d$width, d$height)
      book$pages <- index
      saveRDS(book, paste0(base, ".rds"))
    }
    return(invisible(files))
  }
  d <- draw(book)
  .save_pages(
    d$pages,
    file.path(out_dir, book$folder, paste0(stem, ".pdf")),
    d$width,
    d$height
  )
}

render_figure_books <- function(plot_data, out_dir) {
  for (stem in names(plot_data)) {
    render_figure_from_data(stem, plot_data, out_dir)
  }
  invisible(out_dir)
}
