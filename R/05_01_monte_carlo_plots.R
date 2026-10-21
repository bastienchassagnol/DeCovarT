# Monte Carlo visualisation of simulation benchmarks: long table of
# estimates, raincloud of estimates or errors, forest plot and faceted dots
# of the ADEMP summaries. Algorithm-similarity heatmaps live in
# 05_02_algorithm_similarity.R. Figure books for the article are drawn by
# the scripts under scripts/visualisation/.

#' Pivot Monte Carlo proportion estimates to a long table
#'
#' Turns the `optimisation` block of [run_simulation_benchmark()] into
#' one row per replicate, algorithm, and cell type, with the matching
#' true proportion from `theta_true`. The raincloud intervals drawn from
#' this table are **empirical Monte Carlo quantiles** of
#' \eqn{\hat{\boldsymbol{p}}}, not confidence intervals for
#' \eqn{\boldsymbol{p}} \insertCite{allenRaincloudPlotsMultiplatform2019}{DeCovarT}.
#'
#' @param benchmark List returned by [run_simulation_benchmark()].
#'
#' @return A tibble with `algorithm`, `cell_type`, `estimate`, `p_true`,
#'   `error`, `sample_id`, and any scenario metadata columns.
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
#'     ID = "B1",
#'     true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
#'   ),
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   n = 2,
#'   cores = 1
#' )
#' nrow(pivot_mc_estimates(out))
#' @export
#' @importFrom Rdpack reprompt
#' @seealso [plot_mc_raincloud()], [plot_mc_forest()]
#' @references
#' \insertAllCited{}
#' @srrstats {G2.0} Input is a named list from the benchmark wrapper.
pivot_mc_estimates <- function(benchmark) {
  if (!is.list(benchmark) || is.null(benchmark$optimisation)) {
    stop(
      "`benchmark` must be a list from run_simulation_benchmark().",
      call. = FALSE
    )
  }
  opt <- benchmark$optimisation
  config <- benchmark$config
  truths <- benchmark$theta_true
  if (nrow(opt) == 0L) {
    stop("`benchmark$optimisation` is empty.", call. = FALSE)
  }
  if (length(truths) == 0L) {
    stop("`benchmark$theta_true` is missing.", call. = FALSE)
  }
  if (nrow(config) != length(truths)) {
    stop(
      "`config` and `theta_true` must have one entry per scenario.",
      call. = FALSE
    )
  }
  truth_names <- lapply(truths, function(th) {
    nms <- names(th$p)
    if (is.null(nms)) {
      nms <- colnames(th$mu)
    }
    if (is.null(nms)) {
      nms <- paste0("ct", seq_along(th$p))
    }
    nms
  })
  cell_names <- unique(unlist(truth_names))
  missing_ct <- setdiff(cell_names, names(opt))
  if (length(missing_ct) > 0L) {
    stop(
      "`optimisation` is missing cell-type columns: ",
      toString(missing_ct),
      ".",
      call. = FALSE
    )
  }
  # Scenario keys: atomic columns shared by the design grid and the metrics.
  keys <- setdiff(
    intersect(names(config), names(opt)),
    c("true_parameters", "nobservations")
  )
  keys <- keys[
    !vapply(
      keys,
      function(nm) is.list(config[[nm]]) || is.list(opt[[nm]]),
      logical(1)
    )
  ]
  truth_tbl <- dplyr::bind_rows(lapply(seq_along(truths), function(i) {
    row <- tibble::tibble(
      cell_type = truth_names[[i]],
      p_true = as.numeric(truths[[i]]$p)
    )
    if (length(keys) > 0L) {
      row <- dplyr::bind_cols(config[i, keys, drop = FALSE], row)
    }
    row
  }))
  long <- tidyr::pivot_longer(
    opt,
    cols = dplyr::all_of(cell_names),
    names_to = "cell_type",
    values_to = "estimate"
  )
  long <- dplyr::left_join(long, truth_tbl, by = c(keys, "cell_type"))
  long$error <- long$estimate - long$p_true
  long$algorithm <- .relevel_algorithm(long$algorithm)
  long
}

#' Solver factor with the convolution solvers last
#'
#' Known solvers come first in a fixed display order (Marquardt--Levenberg
#' before Newton--Raphson); other labels follow in order of appearance.
#'
#' @noRd
.relevel_algorithm <- function(x) {
  x <- as.character(x)
  known <- c(
    "nnls",
    "lsei",
    "cibersort",
    "SA",
    "gradient",
    "LBFGS",
    "LBFGSB",
    "L-BFGS-B",
    "Marquardt-Levenberg",
    "Newton-Raphson"
  )
  present <- intersect(known, unique(x))
  factor(x, levels = c(present, setdiff(unique(x), present)))
}

#' Add true-ratio labels (and optional reference lines) per cell type
#'
#' One coloured label per cell type and facet: the most abundant type to the
#' right of its true proportion, the others stacked to the left. With
#' `vlines = TRUE`, a dashed vertical line marks each true proportion.
#' `pal` defaults to four qualitative colours named by cell type.
#'
#' @noRd
.add_true_ratio_labels <- function(
  p,
  df,
  extra_keys = NULL,
  pal = NULL,
  size = 2.6,
  vlines = FALSE
) {
  if (!all(c("cell_type", "p_true") %in% names(df))) {
    return(p)
  }
  keys <- intersect(unique(c("cell_type", extra_keys)), names(df))
  truth <- dplyr::distinct(df, dplyr::across(dplyr::all_of(c(keys, "p_true"))))
  grp <- setdiff(keys, "cell_type")
  place <- function(part) {
    part <- part[order(as.character(part$cell_type)), , drop = FALSE]
    value <- as.numeric(part$p_true)
    abundant <- which.max(value)
    if (length(abundant) != 1L || !is.finite(value[[abundant]])) {
      abundant <- 1L
    }
    right <- seq_len(nrow(part)) == abundant
    rank <- ifelse(right, cumsum(right), cumsum(!right))
    part$lab_x <- pmin(
      pmax(ifelse(right, value + 0.04, value - 0.04), 0.02),
      0.98
    )
    part$lab_y <- Inf
    part$hjust <- ifelse(right, 0, 1)
    part$vjust <- 1.08 + 1.2 * (rank - 1)
    part$lab <- paste0(
      as.character(part$cell_type),
      ": ",
      sprintf("%.2f", value)
    )
    part
  }
  lab_df <- if (length(grp) == 0L) {
    place(truth)
  } else {
    dplyr::bind_rows(lapply(
      split(
        truth,
        interaction(truth[grp], drop = TRUE, lex.order = TRUE),
        drop = TRUE
      ),
      place
    ))
  }
  cell_types <- unique(as.character(lab_df$cell_type))
  if (is.null(pal)) {
    pal <- stats::setNames(
      c("#E41A1C", "#4DAF4A", "#377EB8", "#984EA3")[seq_along(cell_types)],
      cell_types
    )
  }
  for (ct in cell_types) {
    sub <- lab_df[as.character(lab_df$cell_type) == ct, , drop = FALSE]
    col <- unname(pal[ct])
    if (is.na(col) || !nzchar(col)) {
      col <- "#222222"
    }
    if (isTRUE(vlines)) {
      p <- p +
        ggplot2::geom_vline(
          xintercept = unique(sub$p_true),
          colour = col,
          linetype = "longdash",
          linewidth = 0.7
        )
    }
    p <- p +
      ggplot2::geom_label(
        data = sub,
        mapping = ggplot2::aes(
          x = .data[["lab_x"]],
          y = .data[["lab_y"]],
          label = .data[["lab"]],
          hjust = .data[["hjust"]],
          vjust = .data[["vjust"]]
        ),
        inherit.aes = FALSE,
        colour = col,
        fill = ggplot2::alpha("white", 0.92),
        size = size,
        fontface = "bold",
        linewidth = 0.25,
        label.padding = grid::unit(0.12, "lines")
      )
  }
  p
}

#' Optional `facet_grid` from row and column variable names
#'
#' @noRd
.facet_grid_from_names <- function(
  data,
  facet_rows,
  facet_cols,
  scales = "fixed"
) {
  vars <- unique(c(facet_rows, facet_cols))
  vars <- vars[!is.na(vars) & nzchar(vars) & vars != "."]
  unknown <- setdiff(vars, names(data))
  if (length(unknown) > 0L) {
    stop(
      "Facet column '",
      unknown[[1L]],
      "' is not in the plotting table.",
      call. = FALSE
    )
  }
  if (is.null(facet_rows) && is.null(facet_cols)) {
    return(NULL)
  }
  side <- function(v) if (length(v) == 0L) "." else paste(v, collapse = " + ")
  ggplot2::facet_grid(
    stats::as.formula(paste(side(facet_rows), "~", side(facet_cols))),
    labeller = ggplot2::label_value,
    scales = scales,
    drop = FALSE
  )
}

#' Faceted ggplot2 theme (black strips, panel border)
#'
#' `theme_minimal()` plus a panel border and white-on-black facet
#' strips, after `theme_features()` in the
#' [atlas-feature-selection-benchmark](https://github.com/theislab/atlas-feature-selection-benchmark/blob/b89fc0f66747062e6e1b4b35bd392b27ad035295/analysis/R/plotting.R)
#' plotting helpers.
#'
#' @param base_size Base font size for [ggplot2::theme_minimal()].
#' @param ... Passed to [ggplot2::theme()].
#'
#' @return A `ggplot2` theme object.
#' @export
#' @seealso [plot_mc_forest()], [plot_mc_raincloud()]
#' @examples
#' ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) +
#'   ggplot2::geom_point() +
#'   ggplot2::facet_wrap(~cyl) +
#'   theme_decovart_facets()
theme_decovart_facets <- function(base_size = 11, ...) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(
        colour = "grey35",
        fill = NA
      ),
      panel.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      plot.background = ggplot2::element_rect(
        fill = "transparent",
        colour = NA
      ),
      strip.text = ggplot2::element_text(colour = "white"),
      strip.background = ggplot2::element_rect(
        fill = "black",
        colour = NA
      ),
      strip.clip = "off",
      panel.spacing = grid::unit(0.25, "lines"),
      plot.margin = ggplot2::margin(4, 6, 4, 4),
      ...
    )
}

#' Horizontal raincloud of Monte Carlo proportion estimates
#'
#' Half-eye densities, dots, and empirical 50% / 95% intervals
#' \insertCite{allenRaincloudPlotsMultiplatform2019}{DeCovarT} for the
#' Monte Carlo sampling distribution of \eqn{\hat p_j}. Default
#' aesthetics: numeric axis is the estimation error
#' \eqn{\hat p_j-p_j^{\star}} (so a vertical line at 0 is bias);
#' the y-axis is the cell type; colour and dodge group the deconvolution
#' algorithm. Optional `facet_rows` / `facet_cols` split scenarios
#' (for example number of genes versus pairwise cosine).
#'
#' When `quantity = "estimate"` and `slab = "halfeye"`, slabs use a
#' bounded kernel on \eqn{[0,1]} ([ggdist::density_bounded()]) with a
#' coloured outline (`slab_colour` / `slab_linewidth`). Vertex pile-up
#' can leave that KDE empty; `slab = "dotsinterval"`
#' ([ggdist::stat_dotsinterval()]) draws quantile dots instead.
#' True proportions are labelled at the top of each panel (the most
#' abundant type to the right of its reference line; the others to the
#' left).
#'
#' The inner interval is the central 50% of Monte Carlo replicates; the
#' outer interval is the central 95%. These are **not** confidence
#' intervals for \eqn{p_j}.
#'
#' @param benchmark List from [run_simulation_benchmark()], or a table
#'   from [pivot_mc_estimates()].
#' @param quantity `"error"` (default) or `"estimate"`.
#' @param facet_rows,facet_cols Optional column names for
#'   [ggplot2::facet_grid()] rows and columns.
#' @param .width Passed to [ggdist::stat_halfeye()]; default
#'   `c(0.5, 0.95)`.
#' @param include_dots If `FALSE`, omit [ggdist::stat_dots()] (use this
#'   on large Monte Carlo tables: the Fortran QP in `nudge_bins` cannot
#'   take long vectors).
#' @param max_dots Maximum rows passed to `stat_dots` when
#'   `include_dots = TRUE`. Extra rows are subsampled.
#' @param max_rows Optional cap on the plotting table (half-eye and
#'   dots). When the Monte Carlo stack is huge, subsample before
#'   drawing.
#' @param dodge_width Width passed to [ggplot2::position_dodge()].
#' @param slab_scale Slab height for [ggdist::stat_halfeye()].
#' @param slab_alpha Transparency of the density slab (`1` is opaque).
#' @param category_spacing Vertical gap between discrete y categories
#'   (cell types). Values greater than 1 insert extra space so slabs
#'   from neighbouring types do not overlap; the axis still shows the
#'   original labels.
#' @param slab `"halfeye"` (default KDE raincloud) or `"dotsinterval"`
#'   (quantile dots plus interval; avoids empty bounded slabs).
#'
#' @srrstats {G2.3} Restricted character input (`quantity`).
#' @srrstats {G2.3a} Validated via `.match_arg_case_insensitive()`.
#' @srrstats {G2.3b} Matching is case-insensitive.
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
#'     n_genes = 2L,
#'     cosine = 0,
#'     true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
#'   ),
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   n = 4,
#'   cores = 1
#' )
#' if (requireNamespace("ggdist", quietly = TRUE)) {
#'   plot_mc_raincloud(out, facet_cols = "cosine")
#' }
#' @export
#' @seealso [plot_mc_forest()], [pivot_mc_estimates()]
#' @references
#' \insertAllCited{}
plot_mc_raincloud <- function(
  benchmark,
  quantity = c("error", "estimate"),
  facet_rows = NULL,
  facet_cols = NULL,
  .width = c(0.5, 0.95),
  include_dots = TRUE,
  max_dots = 2000L,
  max_rows = NULL,
  dodge_width = 0.95,
  slab_scale = 1.4,
  slab_alpha = 1,
  category_spacing = 1,
  slab = c("halfeye", "dotsinterval")
) {
  .check_suggested_package("ggdist", "plot_mc_raincloud")
  quantity <- .match_arg_case_insensitive(quantity, c("error", "estimate"))
  slab <- .match_arg_case_insensitive(slab, c("halfeye", "dotsinterval"))
  df <- if (is.data.frame(benchmark)) {
    benchmark
  } else {
    pivot_mc_estimates(benchmark)
  }
  if (!is.null(max_rows) && nrow(df) > as.integer(max_rows)) {
    df <- df[sample.int(nrow(df), as.integer(max_rows)), , drop = FALSE]
  }
  needed <- c("algorithm", "cell_type", "estimate", "p_true", "error")
  missing <- setdiff(needed, names(df))
  if (length(missing) > 0L) {
    stop(
      "`benchmark` must contain ",
      toString(missing),
      ".",
      call. = FALSE
    )
  }
  df$cell_type <- factor(df$cell_type, levels = unique(df$cell_type))
  df$algorithm <- .relevel_algorithm(df$algorithm)
  y_levels <- levels(df$cell_type)
  use_spaced_y <- is.finite(category_spacing) &&
    category_spacing > 1 &&
    length(y_levels) > 0L
  if (isTRUE(use_spaced_y)) {
    df$y_plot <- (as.numeric(df$cell_type) - 1) * category_spacing + 1
    y_breaks <- (seq_along(y_levels) - 1) * category_spacing + 1
    y_aes <- "y_plot"
  } else {
    df$y_plot <- df$cell_type
    y_aes <- "cell_type"
  }
  x_var <- if (identical(quantity, "error")) "error" else "estimate"
  x_lab <- if (identical(quantity, "error")) {
    "Monte Carlo error (estimate minus truth)"
  } else {
    "Monte Carlo estimate"
  }
  if (
    identical(quantity, "estimate") &&
      identical(slab, "halfeye") &&
      !isTRUE(include_dots)
  ) {
    grp_cols <- c("algorithm", "cell_type")
    if ("panel" %in% names(df)) {
      grp_cols <- c(grp_cols, "panel")
    }
    nuniq <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(grp_cols))) |>
      dplyr::summarise(
        n_unique = dplyr::n_distinct(round(.data[[x_var]], 8)),
        .groups = "drop"
      )
    if (any(nuniq$n_unique <= 2L, na.rm = TRUE)) {
      include_dots <- TRUE
    }
  }
  if (identical(slab, "dotsinterval")) {
    include_dots <- FALSE
  }
  dodge <- ggplot2::position_dodge(width = dodge_width)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_var]],
      y = .data[[y_aes]],
      fill = .data[["algorithm"]],
      colour = .data[["algorithm"]]
    )
  )
  if (identical(slab, "dotsinterval")) {
    p <- p +
      ggdist::stat_dotsinterval(
        orientation = "horizontal",
        .width = .width,
        quantiles = 100,
        layout = "weave",
        overflow = "compress",
        justification = -0.08,
        point_interval = ggdist::median_qi,
        normalize = "groups",
        scale = slab_scale,
        interval_size = 2.8,
        point_size = 1.6,
        slab_linewidth = 0.4,
        slab_alpha = slab_alpha,
        position = dodge
      )
  } else {
    slab_density <- if (identical(quantity, "estimate")) {
      ggdist::density_bounded(bounds = c(0, 1))
    } else {
      ggdist::density_unbounded()
    }
    slab_limits <- if (identical(quantity, "estimate")) {
      c(0, 1)
    } else {
      NULL
    }
    p <- p +
      ggdist::stat_halfeye(
        ggplot2::aes(slab_colour = ggplot2::after_scale(.data[["fill"]])),
        orientation = "horizontal",
        .width = .width,
        justification = -0.08,
        point_interval = ggdist::median_qi,
        normalize = "groups",
        scale = slab_scale,
        density = slab_density,
        limits = slab_limits,
        interval_size = 2.8,
        point_size = 1.6,
        slab_linewidth = 0.6,
        slab_alpha = slab_alpha,
        position = dodge
      )
  }
  if (isTRUE(include_dots)) {
    df_dots <- df
    n_dots <- nrow(df_dots)
    if (is.finite(max_dots) && n_dots > as.integer(max_dots)) {
      df_dots <- df_dots[
        sample.int(n_dots, as.integer(max_dots)),
        ,
        drop = FALSE
      ]
    }
    p <- p +
      ggdist::stat_dots(
        data = df_dots,
        orientation = "horizontal",
        side = "bottom",
        position = dodge
      )
  }
  p <- p +
    ggplot2::labs(
      x = x_lab,
      y = "Cell type",
      fill = "Algorithm",
      colour = "Algorithm",
      caption = paste(
        "Central 50% and 95% of Monte Carlo replicates;",
        "not a confidence interval for p."
      )
    ) +
    theme_decovart_facets()
  if (isTRUE(use_spaced_y)) {
    p <- p +
      ggplot2::scale_y_continuous(
        breaks = y_breaks,
        labels = y_levels,
        expand = ggplot2::expansion(add = 0.45 * category_spacing)
      )
  }
  if (identical(quantity, "error")) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
  } else {
    extra <- unique(c(facet_rows, facet_cols, "panel", "page"))
    extra <- extra[extra %in% names(df)]
    p <- .add_true_ratio_labels(p, df, extra_keys = extra, vlines = TRUE)
    p <- p +
      ggplot2::coord_cartesian(xlim = c(-0.04, 1.04), clip = "off") +
      ggplot2::labs(
        caption = paste(
          "Central 50% and 95% of Monte Carlo replicates;",
          "not a confidence interval for p.",
          "Dashed vertical lines and labels: true cell-type proportions",
          "(abundant type to the right, others to the left).",
          if (identical(slab, "dotsinterval")) {
            "Dots show 100 quantiles per group (no bounded KDE slab)."
          } else {
            paste(
              "Slabs are bounded to [0, 1] (unit simplex).",
              "When a panel piles on a vertex, dots replace an empty KDE slab."
            )
          }
        )
      )
  }
  facet <- .facet_grid_from_names(df, facet_rows, facet_cols)
  if (!is.null(facet)) {
    p <- p + facet
  }
  p
}

#' Monte Carlo summaries joined with the mean absolute error
#'
#' The `monte_carlo` block of a benchmark plus the MAE of the estimates over
#' replicates, grouped by algorithm, cell type, scenario keys and any extra
#' (facet) columns. Shared by [plot_mc_forest()] and [plot_mc_metric_dots()].
#'
#' @return A list with `table` and the grouping `keys`.
#' @noRd
.mc_summary_table <- function(benchmark, extra_keys = NULL) {
  if (is.data.frame(benchmark)) {
    stop(
      "This plot needs the full run_simulation_benchmark() list.",
      call. = FALSE
    )
  }
  mc <- benchmark$monte_carlo
  long <- pivot_mc_estimates(benchmark)
  keys <- intersect(
    unique(c(intersect(names(mc), names(long)), extra_keys)),
    names(long)
  )
  keys <- keys[!vapply(long[keys], is.list, logical(1))]
  mae <- long |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      mae = mean(abs(.data[["error"]]), na.rm = TRUE),
      .groups = "drop"
    )
  list(
    table = dplyr::left_join(mc, mae, by = intersect(names(mc), keys)),
    keys = keys
  )
}

#' Forest plot of ADEMP Monte Carlo summaries
#'
#' Dot-and-whisker display of bias, RMSE, MAE, coverage, mean interval
#' width, SE/SD ratio, and optimiser failure rate by algorithm and cell
#' type
#' \insertCite{allenRaincloudPlotsMultiplatform2019}{DeCovarT}.
#' Coverage whiskers are the Wilson interval already stored on
#' `monte_carlo` ([coverage_mc_interval()];
#' \insertCite{wilsonProbableInferenceLaw1927}{DeCovarT}): they are
#' intervals for the coverage *rate*, not for \eqn{p_j}. Bias is
#' referenced at 0; coverage at 0.95; SE/SD (`se_sd_ratio`, mean model
#' SE / empirical SD) at 1. Interior Wald SEs come from
#' [vcov_ilr_delta()]. Pairwise algorithm contrasts
#' (MAE differences versus a reference solver on the same Monte Carlo
#' replicates) can be read from the raincloud of paired errors; they do
#' not need a second bootstrap.
#'
#' @inheritParams plot_mc_raincloud
#' @param metrics Character vector of summaries to display.
#'
#' @srrstats {G2.3} Restricted character input (`metrics`).
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
#'     n_genes = 2L,
#'     cosine = 0,
#'     true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
#'   ),
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   n = 4,
#'   cores = 1
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_mc_forest(out, facet_cols = "cosine")
#' }
#' @export
#' @seealso [plot_mc_raincloud()], [coverage_mc_interval()],
#'   [vcov_ilr_delta()], [theme_decovart_facets()]
#' @references
#' \insertAllCited{}
plot_mc_forest <- function(
  benchmark,
  facet_rows = NULL,
  facet_cols = NULL,
  metrics = c(
    "bias",
    "rmse",
    "mae",
    "coverage",
    "mean_interval_width",
    "se_sd_ratio",
    "failure_rate"
  )
) {
  # One row per summary: value column, optional interval columns, reference.
  spec <- tibble::tribble(
    ~metric               , ~value                , ~lo              , ~hi              , ~reference ,
    "bias"                , "bias"                , NA               , NA               , 0          ,
    "rmse"                , "rmse"                , NA               , NA               , NA         ,
    "mae"                 , "mae"                 , NA               , NA               , NA         ,
    "coverage"            , "coverage"            , "coverage_lower" , "coverage_upper" , 0.95       ,
    "mean_interval_width" , "mean_interval_width" , NA               , NA               , NA         ,
    "se_sd_ratio"         , "se_sd_ratio"         , NA               , NA               , 1          ,
    "failure_rate"        , "failure_rate"        , "failure_lower"  , "failure_upper"  , 0
  )
  metrics <- vapply(
    metrics,
    .match_arg_case_insensitive,
    character(1),
    choices = spec$metric
  )
  summary <- .mc_summary_table(benchmark, c(facet_rows, facet_cols))
  forest <- summary$table
  # Optimiser failure rate per algorithm and scenario, with a Wilson interval.
  fail_keys <- intersect(
    unique(c("algorithm", setdiff(summary$keys, "cell_type"))),
    names(benchmark$optimisation)
  )
  fail_tbl <- benchmark$optimisation |>
    dplyr::group_by(dplyr::across(dplyr::all_of(fail_keys))) |>
    dplyr::summarise(
      n_fail = sum(!.data[["numerical_converged"]], na.rm = TRUE),
      n = sum(!is.na(.data[["numerical_converged"]])),
      .groups = "drop"
    )
  fail_int <- purrr::map2(fail_tbl$n_fail, fail_tbl$n, function(s, n) {
    coverage_mc_interval(rep(c(TRUE, FALSE), c(s, n - s)))
  })
  fail_tbl$failure_rate <- vapply(fail_int, `[[`, numeric(1), "coverage")
  fail_tbl$failure_lower <- vapply(fail_int, `[[`, numeric(1), "lower")
  fail_tbl$failure_upper <- vapply(fail_int, `[[`, numeric(1), "upper")
  forest <- dplyr::left_join(
    forest,
    fail_tbl[, setdiff(names(fail_tbl), c("n_fail", "n"))],
    by = intersect(fail_keys, names(forest))
  )
  pieces <- lapply(intersect(spec$metric, metrics), function(m) {
    s <- spec[spec$metric == m, ]
    rows <- forest
    if (identical(m, "failure_rate")) {
      keep <- intersect(
        c(fail_keys, "failure_rate", "failure_lower", "failure_upper"),
        names(forest)
      )
      rows <- dplyr::distinct(forest, dplyr::across(dplyr::all_of(keep)))
      rows$cell_type <- "all"
    }
    tibble::tibble(
      rows,
      metric = m,
      estimate = rows[[s$value]],
      lower = if (is.na(s$lo)) NA_real_ else rows[[s$lo]],
      upper = if (is.na(s$hi)) NA_real_ else rows[[s$hi]],
      reference = s$reference
    )
  })
  keep <- unique(c(
    "algorithm",
    "cell_type",
    "metric",
    "estimate",
    "lower",
    "upper",
    "reference",
    facet_rows,
    facet_cols
  ))
  plot_df <- dplyr::bind_rows(lapply(pieces, function(x) {
    dplyr::select(x, dplyr::any_of(keep))
  }))
  plot_df <- dplyr::filter(plot_df, is.finite(.data[["estimate"]]))
  plot_df$metric <- factor(plot_df$metric, levels = unique(metrics))
  vline_df <- dplyr::distinct(
    dplyr::filter(plot_df, is.finite(.data[["reference"]])),
    .data[["metric"]],
    .data[["reference"]]
  )
  interval_df <- dplyr::filter(
    plot_df,
    is.finite(.data[["lower"]]) & is.finite(.data[["upper"]])
  )
  dodge <- ggplot2::position_dodge(width = 0.6)
  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["estimate"]],
      y = .data[["algorithm"]],
      colour = .data[["cell_type"]]
    )
  ) +
    ggplot2::geom_vline(
      data = vline_df,
      ggplot2::aes(xintercept = .data[["reference"]]),
      linetype = "dashed",
      colour = "grey40"
    ) +
    ggplot2::geom_linerange(
      data = interval_df,
      ggplot2::aes(xmin = .data[["lower"]], xmax = .data[["upper"]]),
      position = dodge
    ) +
    ggplot2::geom_point(position = dodge) +
    ggplot2::labs(
      x = "Monte Carlo summary",
      y = "Algorithm",
      colour = "Cell type",
      caption = paste(
        "Coverage whiskers are Wilson intervals for the coverage rate;",
        "bias reference is 0; coverage reference is 0.95;",
        "SE/SD is mean model SE / empirical SD (1 = calibrated).",
        "Interior Wald SEs use the ILR expected-Fisher delta method."
      )
    ) +
    theme_decovart_facets()
  facet <- .facet_grid_from_names(
    plot_df,
    facet_rows,
    unique(c("metric", facet_cols)),
    scales = "free_x"
  )
  if (!is.null(facet)) {
    p <- p + facet
  }
  p
}

#' Faceted dot plot of several ADEMP metrics
#'
#' One panel per metric; colour is a **min-max normalised** score with
#' \(1\) = best inside that metric (and facet). Do not map a second
#' primary score to point size. Coverage is scored as
#' \eqn{|\widehat C-0.95|}, so over-wide intervals that inflate
#' coverage are not rewarded. There is no default weighted composite.
#'
#' @inheritParams plot_mc_forest
#' @param weights Optional named numeric weights (sum to 1) used only
#'   to add a `composite` facet beside the components.
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
#'     cosine = 0,
#'     true_theta = list(list(p = c(0.5, 0.5), mu = mu, sigma = Sigma))
#'   ),
#'   deconvolution_functions = list(
#'     "nnls" = list(FUN = deconvolute_ratios_nnls)
#'   ),
#'   n = 4,
#'   cores = 1
#' )
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_mc_metric_dots(out, facet_cols = "cosine")
#' }
#' @export
#' @seealso [plot_mc_forest()], [plot_algorithm_similarity()]
#' @srrstats {G2.3} Restricted character input (`metrics`).
plot_mc_metric_dots <- function(
  benchmark,
  facet_rows = NULL,
  facet_cols = NULL,
  metrics = c("rmse", "mae", "coverage", "mean_interval_width"),
  weights = NULL
) {
  # Raw score of each metric, oriented so that lower is better.
  raw_of <- list(
    bias = function(x) abs(x$bias),
    rmse = function(x) x$rmse,
    mae = function(x) x$mae,
    coverage = function(x) abs(x$coverage - 0.95),
    mean_interval_width = function(x) x$mean_interval_width
  )
  metrics <- vapply(
    metrics,
    .match_arg_case_insensitive,
    character(1),
    choices = names(raw_of)
  )
  summary <- .mc_summary_table(benchmark)
  plot_df <- dplyr::bind_rows(lapply(
    intersect(names(raw_of), metrics),
    function(m) {
      rows <- summary$table
      dplyr::select(
        tibble::tibble(rows, metric = m, raw = raw_of[[m]](rows)),
        dplyr::any_of(c(summary$keys, "metric", "raw"))
      )
    }
  ))
  # Min-max rescaling inside each metric and facet: 1 = best.
  rescale <- function(x) {
    finite <- x[is.finite(x)]
    if (length(finite) == 0L) {
      return(rep(NA_real_, length(x)))
    }
    rng <- max(finite) - min(finite)
    if (rng < .Machine$double.eps) {
      return(ifelse(is.finite(x), 1, NA_real_))
    }
    1 - (x - min(finite)) / rng
  }
  scale_cols <- intersect(
    unique(c("metric", facet_rows, facet_cols)),
    names(plot_df)
  )
  plot_df <- plot_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scale_cols))) |>
    dplyr::mutate(score = rescale(.data[["raw"]])) |>
    dplyr::ungroup()
  if (!is.null(weights)) {
    w_names <- names(weights)
    if (is.null(w_names) || any(w_names == "")) {
      stop("`weights` must be a named numeric vector.", call. = FALSE)
    }
    weights <- weights / sum(weights)
    missing_w <- setdiff(w_names, unique(as.character(plot_df$metric)))
    if (length(missing_w) > 0L) {
      stop("Unknown weight names: ", toString(missing_w), ".", call. = FALSE)
    }
    keys <- setdiff(names(plot_df), c("metric", "raw", "score"))
    composite <- plot_df |>
      dplyr::filter(.data[["metric"]] %in% w_names) |>
      dplyr::mutate(w = unname(weights[as.character(.data[["metric"]])])) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(
        metric = "composite",
        raw = NA_real_,
        score = sum(.data[["w"]] * .data[["score"]], na.rm = TRUE),
        .groups = "drop"
      )
    plot_df <- dplyr::bind_rows(plot_df, composite)
  }
  plot_df$metric <- factor(
    plot_df$metric,
    levels = unique(c(metrics, if (!is.null(weights)) "composite"))
  )
  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = .data[["cell_type"]],
      y = .data[["algorithm"]],
      colour = .data[["score"]]
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_colour_gradient(
      low = "#d73027",
      high = "#1a9850",
      limits = c(0, 1),
      name = "Relative score (1 = best)"
    ) +
    ggplot2::labs(
      x = "Cell type",
      y = "Algorithm",
      caption = paste(
        "Colour is min-max scaled within each metric;",
        "coverage uses |C - 0.95|. No default composite score."
      )
    ) +
    theme_decovart_facets() +
    .facet_grid_from_names(plot_df, facet_rows, unique(c("metric", facet_cols)))
}
