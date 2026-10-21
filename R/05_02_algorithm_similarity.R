# Behavioural similarity of deconvolution algorithms: Pearson correlation of
# their Monte Carlo estimates and a Ward-ordered tile heatmap with an
# optional dendrogram.

#' Algorithm-similarity correlation from a Monte Carlo benchmark
#'
#' Pearson correlations
#' \eqn{r_{ab}=\mathrm{cor}(\hat p_a,\hat p_b)} across Monte Carlo
#' replicates (and cell types). This is **behavioural similarity**: two
#' solvers can correlate near 1 while remaining systematically biased.
#' Hierarchical order uses Ward D2 on \eqn{d_{ab} = 1 - r_{ab}}.
#'
#' @param benchmark List from [run_simulation_benchmark()].
#' @param facet_rows,facet_cols Optional scenario columns. When supplied,
#'   correlations are computed inside each facet cell.
#'
#' @return A tibble with `algorithm_x`, `algorithm_y`, `correlation`,
#'   and any facet columns.
#'
#' @examplesIf .Platform$OS.type != "windows"
#' set.seed(1)
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' out <- run_simulation_benchmark(
#'   tibble::tibble(
#'     true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
#'   ),
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls),
#'     "rlm" = list(FUN = deconvolute_ratios_rlm)
#'   ),
#'   n = 4,
#'   cores = 1
#' )
#' algorithm_similarity(out)
#' @export
#' @seealso [plot_algorithm_similarity()], [pivot_mc_estimates()]
#' @srrstats {G2.0} Requires a named benchmark list with at least two
#'   algorithms.
algorithm_similarity <- function(
  benchmark,
  facet_rows = NULL,
  facet_cols = NULL
) {
  long <- pivot_mc_estimates(benchmark)
  algos <- levels(droplevels(long$algorithm))
  if (length(algos) < 2L) {
    stop(
      "algorithm_similarity() needs at least two algorithms.",
      call. = FALSE
    )
  }
  group_cols <- intersect(unique(c(facet_rows, facet_cols)), names(long))
  id_cols <- intersect(c("ID", "sample_id", "cell_type"), names(long))
  groups <- if (length(group_cols) == 0L) {
    list(long)
  } else {
    split(long, interaction(long[, group_cols, drop = FALSE], drop = TRUE))
  }
  dplyr::bind_rows(lapply(groups, function(part) {
    wide <- tidyr::pivot_wider(
      part,
      id_cols = dplyr::all_of(id_cols),
      names_from = "algorithm",
      values_from = "estimate"
    )
    present <- intersect(algos, names(wide))
    if (length(present) < 2L) {
      stop(
        "Need at least two algorithms with paired Monte Carlo estimates.",
        call. = FALSE
      )
    }
    r_mat <- stats::cor(
      as.matrix(wide[, present, drop = FALSE]),
      use = "pairwise.complete.obs"
    )
    tbl <- tibble::tibble(
      algorithm_x = rep(present, times = length(present)),
      algorithm_y = rep(present, each = length(present)),
      correlation = as.vector(r_mat)
    )
    if (length(group_cols) > 0L) {
      tbl <- dplyr::bind_cols(part[1L, group_cols, drop = FALSE], tbl)
    }
    tbl
  }))
}

#' Tile heatmap of algorithm-similarity correlations
#'
#' `ggplot2::geom_tile()` display of [algorithm_similarity()], with
#' rows and columns ordered by Ward D2 clustering of
#' \(1-r\). Optional dendrogram via `ggdendro` (Suggests), drawn to
#' the **right** of the tiles with leaves flush against the heatmap
#' and Ward merge heights on the axis.
#'
#' @inheritParams algorithm_similarity
#' @param dendrogram If `TRUE`, attach a `ggdendro` ggplot to the right
#'   of the tiles as attribute `"dendrogram"` (ignored when scenario
#'   facets are used).
#'
#' @return A `ggplot` object.
#'
#' @examplesIf .Platform$OS.type != "windows"
#' set.seed(1)
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' out <- run_simulation_benchmark(
#'   tibble::tibble(
#'     true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
#'   ),
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls),
#'     "rlm" = list(FUN = deconvolute_ratios_rlm)
#'   ),
#'   n = 4,
#'   cores = 1
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_algorithm_similarity(out)
#' }
#' @export
#' @seealso [algorithm_similarity()], [plot_mc_metric_dots()]
plot_algorithm_similarity <- function(
  benchmark,
  facet_rows = NULL,
  facet_cols = NULL,
  dendrogram = FALSE
) {
  sim <- algorithm_similarity(
    benchmark,
    facet_rows = facet_rows,
    facet_cols = facet_cols
  )
  group_cols <- unique(c(facet_rows, facet_cols))
  group_cols <- group_cols[group_cols %in% names(sim)]
  order_src <- sim
  if (length(group_cols) > 0L) {
    order_src <- sim |>
      dplyr::group_by(.data[["algorithm_x"]], .data[["algorithm_y"]]) |>
      dplyr::summarise(
        correlation = mean(.data[["correlation"]], na.rm = TRUE),
        .groups = "drop"
      )
  }
  r_mat <- unclass(stats::xtabs(
    correlation ~ algorithm_x + algorithm_y,
    data = order_src
  ))
  # Ward D2 on the correlation distance 1 - r (clamped at 0).
  hc <- stats::hclust(stats::as.dist(pmax(1 - r_mat, 0)), method = "ward.D2")
  ord <- hc$labels[hc$order]
  sim$algorithm_x <- factor(sim$algorithm_x, levels = ord)
  sim$algorithm_y <- factor(sim$algorithm_y, levels = rev(ord))
  add_dend <- isTRUE(dendrogram) && length(group_cols) == 0L
  p <- ggplot2::ggplot(
    sim,
    ggplot2::aes(
      x = .data[["algorithm_x"]],
      y = .data[["algorithm_y"]],
      fill = .data[["correlation"]]
    )
  ) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_gradient2(
      limits = c(-1, 1),
      midpoint = 0,
      low = "#2166ac",
      mid = "white",
      high = "#b2182b"
    ) +
    ggplot2::labs(
      x = "Algorithm",
      y = "Algorithm",
      fill = "Pearson r"
    ) +
    theme_decovart_facets() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5
      ),
      legend.position = "bottom"
    )
  if (isTRUE(add_dend)) {
    p <- p +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_cartesian(expand = FALSE)
  } else {
    p <- p + ggplot2::coord_equal()
  }
  facet <- .facet_grid_from_names(sim, facet_rows, facet_cols)
  if (!is.null(facet)) {
    p <- p + facet
  }
  if (isTRUE(dendrogram)) {
    if (length(group_cols) > 0L) {
      warning(
        "dendrogram = TRUE is ignored when scenario facets are used.",
        call. = FALSE
      )
    } else {
      .check_suggested_package("ggdendro", "plot_algorithm_similarity")
      .check_suggested_package("cowplot", "plot_algorithm_similarity")
      dend <- .similarity_dendrogram_plot(hc, n_leaf = length(ord))
      tiles <- p +
        ggplot2::theme(
          plot.margin = ggplot2::margin(4, 0, 4, 4),
          legend.position = "none"
        )
      combined <- cowplot::plot_grid(
        tiles,
        dend,
        nrow = 1,
        rel_widths = c(1, 0.28),
        align = "h",
        axis = "tb"
      )
      attr(combined, "dendrogram") <- dend
      return(combined)
    }
  }
  p
}

#' Horizontal dendrogram with leaves flush to the left (heatmap side)
#'
#' @noRd
.similarity_dendrogram_plot <- function(hc, n_leaf) {
  ddata <- ggdendro::dendro_data(
    stats::as.dendrogram(hc),
    type = "rectangle"
  )
  seg <- ggdendro::segment(ddata)
  # Leaf 1 of hclust order sits at the top of the heatmap
  # (`algorithm_y = rev(ord)`). Flip the dendrogram index to match.
  seg$x <- n_leaf + 1 - seg$x
  seg$xend <- n_leaf + 1 - seg$xend
  merge_seg <- dplyr::filter(
    seg,
    abs(.data[["y"]] - .data[["yend"]]) < 1e-8,
    abs(.data[["x"]] - .data[["xend"]]) > 1e-8
  )
  merge_lab <- dplyr::mutate(
    merge_seg,
    lab_x = .data[["y"]],
    lab_y = (.data[["x"]] + .data[["xend"]]) / 2,
    lab = sprintf("%.2f", .data[["y"]])
  )
  merge_lab <- dplyr::distinct(
    merge_lab,
    .data[["lab"]],
    .data[["lab_y"]],
    .keep_all = TRUE
  )
  ggplot2::ggplot(seg) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = .data[["y"]],
        y = .data[["x"]],
        xend = .data[["yend"]],
        yend = .data[["xend"]]
      ),
      linewidth = 0.45,
      colour = "grey20",
      lineend = "square"
    ) +
    ggplot2::geom_text(
      data = merge_lab,
      ggplot2::aes(
        x = .data[["lab_x"]],
        y = .data[["lab_y"]],
        label = .data[["lab"]]
      ),
      size = 2.1,
      colour = "grey15",
      hjust = -0.15,
      vjust = 0.5,
      inherit.aes = FALSE
    ) +
    ggplot2::labs(x = "Ward D2 height", y = NULL) +
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0, 0.18))
    ) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, n_leaf + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 8),
      axis.ticks.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 6
      ),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      plot.margin = ggplot2::margin(4, 8, 4, 0)
    )
}
