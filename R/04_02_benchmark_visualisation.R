#' Check optional heatmap Suggests packages
#'
#' @keywords internal
#' @noRd
.check_heatmap_dependencies <- function() {
  pkgs <- c("ComplexHeatmap", "circlize", "viridis")
  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  if (length(missing)) {
    stop(
      "plot_correlation_Heatmap() requires the optional package",
      if (length(missing) > 1L) "s " else " ",
      toString(paste0("'", missing, "'")),
      ". Install CRAN dependencies with install.packages(), and ",
      "'ComplexHeatmap' with BiocManager::install('ComplexHeatmap').",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Plot deconvolution metric heatmaps
#'
#' @description
#' For each algorithm, visualises a selected score over the design grid of
#' simulated \eqn{\boldsymbol{p}} (and related scenario factors).
#' Requires optional Suggests packages `ComplexHeatmap`, `circlize`, and
#' `viridis` (install `ComplexHeatmap` via Bioconductor). Use this
#' helper only for **linked, annotated, hierarchical** grids. For
#' algorithm-similarity correlation matrices prefer
#' [plot_algorithm_similarity()] (`ggplot2::geom_tile()`).
#'
#' @param distribution_metrics Tibble of metric scores from a benchmark.
#' @param score_variable Column name of the metric to display
#'   (`"model_mse"`, `"model_rmse"`, `"model_coef_determination"`,
#'   `"model_coef_determination_adjusted"`, `"model_mae"`, `"model_cor"`,
#'   `"model_ccc"`). Matching is case-insensitive.
#'
#' @srrstats {G2.3} Restricted character input (`score_variable`).
#' @srrstats {G2.3a} Validated via `.match_arg_case_insensitive()` (a `match.arg()`
#'   equivalent).
#' @srrstats {G2.3b} Matching is case-insensitive (`tolower()`).
#' @param n_break Number of colour breaks.
#' @param uni_scale If `FALSE`, each panel uses its own colour scale.
#' @param file Optional PDF path. When supplied, heatmaps are drawn with
#'   [grDevices::pdf()]; a missing `.pdf` suffix is added (G4.0).
#'
#' @return A named list of `ComplexHeatmap` heatmap objects (one per algorithm).
#'
#' @srrstats {G4.0} `file` is passed through `.ensure_file_suffix()`.
#'
#' @examples
#' metrics <- tibble::tibble(
#'   correlation_celltype1 = c(0, 0, 0.5, 0.5),
#'   correlation_celltype2 = c(0, 0.5, 0, 0.5),
#'   algorithm = "nnls",
#'   model_mse = c(0.01, 0.02, 0.015, 0.03)
#' )
#' if (
#'   requireNamespace("ComplexHeatmap", quietly = TRUE) &&
#'     requireNamespace("circlize", quietly = TRUE) &&
#'     requireNamespace("viridis", quietly = TRUE)
#' ) {
#'   ht <- plot_correlation_Heatmap(metrics, score_variable = "model_mse")
#'   names(ht)
#' }
#' @importFrom rlang .data
#' @seealso [plot_algorithm_similarity()], [plot_mc_metric_dots()]
#' @export
plot_correlation_Heatmap <- function(
  distribution_metrics,
  score_variable = "model_mse",
  n_break = 20,
  uni_scale = TRUE,
  file = NULL
) {
  .check_heatmap_dependencies()
  score_variable <- .match_arg_case_insensitive(
    score_variable,
    c(
      "model_mse",
      "model_rmse",
      "model_coef_determination",
      "model_coef_determination_adjusted",
      "model_mae",
      "model_cor",
      "model_ccc"
    )
  )
  distribution_metrics <- distribution_metrics |>
    dplyr::select(dplyr::all_of(c(
      "correlation_celltype1",
      "correlation_celltype2",
      "algorithm",
      score_variable
    ))) |>
    dplyr::mutate(
      algorithm = factor(
        .data[["algorithm"]],
        levels = unique(.data[["algorithm"]])
      )
    )

  # design the scaling colour
  mean_distribution_metrics <- distribution_metrics |>
    dplyr::group_by(
      .data[["correlation_celltype1"]],
      .data[["correlation_celltype2"]],
      .data[["algorithm"]]
    ) |>
    dplyr::summarise(
      mean_metric = mean(.data[[score_variable]], na.rm = TRUE)
    )

  metric_vals <- dplyr::pull(mean_distribution_metrics, "mean_metric")
  min_metric <- min(metric_vals, na.rm = TRUE)
  max_metric <- max(metric_vals, na.rm = TRUE)
  if (uni_scale) {
    col <- circlize::colorRamp2(
      seq(min_metric, max_metric, length.out = n_break),
      viridis::viridis(n_break)
    )
  }

  complex_heatmap_list <- purrr::imap(
    split(distribution_metrics, distribution_metrics[["algorithm"]]),
    function(.mean_signature_matrix, .y) {
      cor_matrix_per_algo <- .mean_signature_matrix |>
        dplyr::select(-"algorithm") |>
        tidyr::pivot_wider(
          names_from = dplyr::all_of("correlation_celltype2"),
          values_from = dplyr::all_of(score_variable),
          values_fn = mean
        ) |>
        tibble::column_to_rownames("correlation_celltype1") |>
        as.matrix()

      if (!uni_scale) {
        col <- circlize::colorRamp2(
          c(
            min(cor_matrix_per_algo, na.rm = TRUE),
            stats::median(cor_matrix_per_algo, na.rm = TRUE),
            max(cor_matrix_per_algo, na.rm = TRUE)
          ),
          c("blue", "white", "red")
        )
      }

      score_label <- toupper(gsub(
        "model_",
        "",
        score_variable,
        fixed = TRUE
      ))
      complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(
        cor_matrix_per_algo,
        col = col,
        name = score_label,
        heatmap_legend_param = list(title = score_label),
        row_title = "Corr cell type 1",
        cluster_rows = FALSE,
        row_names_gp = grid::gpar(fontsize = 8),
        row_labels = colnames(cor_matrix_per_algo),
        row_title_gp = grid::gpar(fontsize = 10),
        column_names_rot = 0,
        cluster_columns = FALSE,
        column_names_gp = grid::gpar(fontsize = 8),
        column_title_gp = grid::gpar(fontsize = 10),
        column_title = "Corr cell type 2",
        column_labels = colnames(cor_matrix_per_algo),
        width = grid::unit(6, "cm"),
        height = grid::unit(6, "cm")
      )
      return(complex_heatmap_per_algo)
    }
  )
  if (!is.null(file)) {
    .write_artefact(complex_heatmap_list, file, kind = "pdf")
  }
  return(complex_heatmap_list)
}

#' Check optional ggplot2 / ggdist Suggests packages
#'
#' @keywords internal
#' @noRd
.check_plot_dependencies <- function(
  need_ggdist = FALSE,
  need_ggdendro = FALSE
) {
  pkgs <- "ggplot2"
  if (isTRUE(need_ggdist)) {
    pkgs <- c(pkgs, "ggdist")
  }
  if (isTRUE(need_ggdendro)) {
    pkgs <- c(pkgs, "ggdendro")
  }
  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    )
  ]
  if (length(missing)) {
    stop(
      "This plot requires the optional package",
      if (length(missing) > 1L) "s " else " ",
      toString(paste0("'", missing, "'")),
      ". Install with install.packages().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Atomic scenario keys shared by config and optimisation tables
#'
#' @noRd
.benchmark_meta_keys <- function(config, opt) {
  shared <- intersect(names(config), names(opt))
  shared <- setdiff(shared, c("true_parameters", "nobservations"))
  keep <- vapply(
    shared,
    function(nm) {
      !is.list(config[[nm]]) && !is.list(opt[[nm]])
    },
    logical(1)
  )
  shared[keep]
}

#' Cell-type names from a `theta_true` element
#'
#' @noRd
.truth_cell_names <- function(th) {
  p <- th$p
  nms <- names(p)
  if (!is.null(nms)) {
    return(nms)
  }
  cn <- colnames(th$mu)
  if (!is.null(cn)) {
    return(cn)
  }
  paste0("ct", seq_along(p))
}

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
#' @examples
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
  if (is.null(truths) || length(truths) == 0L) {
    stop("`benchmark$theta_true` is missing.", call. = FALSE)
  }
  cell_names <- unique(unlist(purrr::map(truths, .truth_cell_names)))
  missing_ct <- setdiff(cell_names, names(opt))
  if (length(missing_ct) > 0L) {
    stop(
      "`optimisation` is missing cell-type columns: ",
      toString(missing_ct),
      ".",
      call. = FALSE
    )
  }
  keys <- .benchmark_meta_keys(config, opt)
  if (nrow(config) != length(truths)) {
    stop(
      "`config` and `theta_true` must have one entry per scenario.",
      call. = FALSE
    )
  }
  truth_tbl <- purrr::map(
    seq_along(truths),
    function(i) {
      th <- truths[[i]]
      nms <- .truth_cell_names(th)
      row <- tibble::tibble(
        cell_type = nms,
        p_true = as.numeric(th$p)
      )
      if (length(keys) > 0L) {
        meta <- config[i, keys, drop = FALSE]
        dplyr::bind_cols(meta, row)
      } else {
        row
      }
    }
  )
  truth_tbl <- dplyr::bind_rows(truth_tbl)
  long <- tidyr::pivot_longer(
    opt,
    cols = dplyr::all_of(cell_names),
    names_to = "cell_type",
    values_to = "estimate"
  )
  if (length(keys) > 0L) {
    long <- dplyr::left_join(long, truth_tbl, by = c(keys, "cell_type"))
  } else {
    long <- dplyr::left_join(long, truth_tbl, by = "cell_type")
  }
  long$error <- long$estimate - long$p_true
  long
}

#' Optional `facet_grid` from row and column variable names
#'
#' @noRd
.facet_grid_from_names <- function(data, facet_rows, facet_cols) {
  for (nm in c(facet_rows, facet_cols)) {
    if (!is.null(nm) && !nm %in% names(data)) {
      stop(
        "Facet column '",
        nm,
        "' is not in the plotting table.",
        call. = FALSE
      )
    }
  }
  if (is.null(facet_rows) && is.null(facet_cols)) {
    return(NULL)
  }
  rows <- if (is.null(facet_rows)) "." else facet_rows
  cols <- if (is.null(facet_cols)) "." else facet_cols
  ggplot2::facet_grid(
    stats::as.formula(paste(rows, "~", cols)),
    labeller = ggplot2::label_both
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
#'
#' @srrstats {G2.3} Restricted character input (`quantity`).
#' @srrstats {G2.3a} Validated via `.match_arg_case_insensitive()`.
#' @srrstats {G2.3b} Matching is case-insensitive.
#'
#' @return A `ggplot` object.
#'
#' @examples
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
  .width = c(0.5, 0.95)
) {
  .check_plot_dependencies(need_ggdist = TRUE)
  quantity <- .match_arg_case_insensitive(quantity, c("error", "estimate"))
  df <- if (is.data.frame(benchmark)) {
    benchmark
  } else {
    pivot_mc_estimates(benchmark)
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
  df$algorithm <- factor(df$algorithm, levels = unique(df$algorithm))
  x_var <- if (identical(quantity, "error")) "error" else "estimate"
  x_lab <- if (identical(quantity, "error")) {
    "Monte Carlo error (estimate minus truth)"
  } else {
    "Monte Carlo estimate"
  }
  dodge <- ggplot2::position_dodge(width = 0.75)
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(
      x = .data[[x_var]],
      y = .data[["cell_type"]],
      fill = .data[["algorithm"]],
      colour = .data[["algorithm"]]
    )
  ) +
    ggdist::stat_halfeye(
      orientation = "horizontal",
      .width = .width,
      justification = -0.15,
      point_interval = ggdist::median_qi,
      normalize = "groups",
      position = dodge
    ) +
    ggdist::stat_dots(
      orientation = "horizontal",
      side = "bottom",
      position = dodge
    ) +
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
    ggplot2::theme_bw()
  if (identical(quantity, "error")) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
  } else {
    truth_cols <- unique(c("cell_type", "p_true", facet_rows, facet_cols))
    truth_cols <- truth_cols[truth_cols %in% names(df)]
    truth <- dplyr::distinct(df, dplyr::across(dplyr::all_of(truth_cols)))
    p <- p +
      ggplot2::geom_point(
        data = truth,
        ggplot2::aes(
          x = .data[["p_true"]],
          y = .data[["cell_type"]]
        ),
        inherit.aes = FALSE,
        shape = 124,
        size = 4
      )
  }
  facet <- .facet_grid_from_names(df, facet_rows, facet_cols)
  if (!is.null(facet)) {
    p <- p + facet
  }
  p
}

#' Wilson interval from event counts
#'
#' @noRd
.wilson_from_counts <- function(n_event, n) {
  purrr::map2(
    as.integer(n_event),
    as.integer(n),
    function(s, nn) {
      s <- max(s, 0L)
      nn <- max(nn, 0L)
      s <- min(s, nn)
      coverage_mc_interval(c(rep(TRUE, s), rep(FALSE, nn - s)))
    }
  )
}

#' Forest plot of ADEMP Monte Carlo summaries
#'
#' Dot-and-whisker display of bias, RMSE, MAE, coverage, mean interval
#' width, and optimiser failure rate by algorithm and cell type
#' \insertCite{allenRaincloudPlotsMultiplatform2019}{DeCovarT}.
#' Coverage whiskers are the Wilson interval already stored on
#' `monte_carlo` ([coverage_mc_interval()];
#' \insertCite{wilsonProbableInferenceLaw1927}{DeCovarT}): they are
#' intervals for the coverage *rate*, not for \eqn{p_j}. Bias is
#' referenced at 0; coverage at 0.95. Pairwise algorithm contrasts
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
#' @examples
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
#' @seealso [plot_mc_raincloud()], [coverage_mc_interval()]
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
    "failure_rate"
  )
) {
  .check_plot_dependencies(need_ggdist = FALSE)
  allowed <- c(
    "bias",
    "rmse",
    "mae",
    "coverage",
    "mean_interval_width",
    "failure_rate"
  )
  metrics <- vapply(
    metrics,
    function(m) .match_arg_case_insensitive(m, allowed),
    character(1)
  )
  if (is.data.frame(benchmark)) {
    stop(
      "plot_mc_forest() needs the full run_simulation_benchmark() list.",
      call. = FALSE
    )
  }
  mc <- benchmark$monte_carlo
  long_est <- pivot_mc_estimates(benchmark)
  keys <- unique(c(
    "algorithm",
    "cell_type",
    .benchmark_meta_keys(benchmark$config, benchmark$optimisation),
    facet_rows,
    facet_cols
  ))
  keys <- keys[keys %in% names(long_est)]
  mae_tbl <- long_est |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      mae = mean(abs(.data[["error"]]), na.rm = TRUE),
      .groups = "drop"
    )
  fail_keys <- unique(c("algorithm", setdiff(keys, "cell_type")))
  fail_keys <- fail_keys[fail_keys %in% names(benchmark$optimisation)]
  fail_tbl <- benchmark$optimisation |>
    dplyr::group_by(dplyr::across(dplyr::all_of(fail_keys))) |>
    dplyr::summarise(
      n_fail = sum(!.data[["numerical_converged"]], na.rm = TRUE),
      n = sum(!is.na(.data[["numerical_converged"]])),
      .groups = "drop"
    )
  fail_int <- .wilson_from_counts(fail_tbl$n_fail, fail_tbl$n)
  fail_tbl$failure_rate <- vapply(
    fail_int,
    `[[`,
    numeric(1),
    "coverage"
  )
  fail_tbl$failure_lower <- vapply(fail_int, `[[`, numeric(1), "lower")
  fail_tbl$failure_upper <- vapply(fail_int, `[[`, numeric(1), "upper")
  fail_tbl$n_fail <- NULL
  fail_tbl$n <- NULL
  forest <- dplyr::left_join(mc, mae_tbl, by = intersect(names(mc), keys))
  forest <- dplyr::left_join(
    forest,
    fail_tbl,
    by = intersect(names(fail_tbl), names(forest))
  )
  pieces <- list()
  if ("bias" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "bias",
        estimate = .data[["bias"]],
        lower = NA_real_,
        upper = NA_real_,
        reference = 0
      )
  }
  if ("rmse" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "rmse",
        estimate = .data[["rmse"]],
        lower = NA_real_,
        upper = NA_real_,
        reference = NA_real_
      )
  }
  if ("mae" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "mae",
        estimate = .data[["mae"]],
        lower = NA_real_,
        upper = NA_real_,
        reference = NA_real_
      )
  }
  if ("coverage" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "coverage",
        estimate = .data[["coverage"]],
        lower = .data[["coverage_lower"]],
        upper = .data[["coverage_upper"]],
        reference = 0.95
      )
  }
  if ("mean_interval_width" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "mean_interval_width",
        estimate = .data[["mean_interval_width"]],
        lower = NA_real_,
        upper = NA_real_,
        reference = NA_real_
      )
  }
  if ("failure_rate" %in% metrics) {
    fail_keep <- unique(c(
      fail_keys,
      "failure_rate",
      "failure_lower",
      "failure_upper"
    ))
    fail_keep <- fail_keep[fail_keep %in% names(forest)]
    fail_plot <- dplyr::distinct(
      forest,
      dplyr::across(dplyr::all_of(fail_keep))
    )
    pieces[[length(pieces) + 1L]] <- fail_plot |>
      dplyr::mutate(
        metric = "failure_rate",
        cell_type = "all",
        estimate = .data[["failure_rate"]],
        lower = .data[["failure_lower"]],
        upper = .data[["failure_upper"]],
        reference = 0
      )
  }
  if (length(pieces) == 0L) {
    stop("`metrics` must contain at least one known summary.", call. = FALSE)
  }
  keep <- c(
    "algorithm",
    "cell_type",
    "metric",
    "estimate",
    "lower",
    "upper",
    "reference",
    facet_rows,
    facet_cols
  )
  keep <- unique(keep[!is.null(keep) & keep %in% names(pieces[[1L]])])
  selected <- purrr::map(
    pieces,
    function(x) dplyr::select(x, dplyr::any_of(keep))
  )
  plot_df <- dplyr::bind_rows(selected)
  plot_df$metric <- factor(plot_df$metric, levels = unique(metrics))
  vline_src <- dplyr::filter(
    plot_df,
    is.finite(.data[["reference"]])
  )
  vline_df <- dplyr::distinct(
    vline_src,
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
        "bias reference is 0; coverage reference is 0.95."
      )
    ) +
    ggplot2::theme_bw()
  rows <- if (is.null(facet_rows)) "." else facet_rows
  cols <- if (is.null(facet_cols)) {
    "metric"
  } else {
    paste("metric", "+", facet_cols)
  }
  p <- p +
    ggplot2::facet_grid(
      stats::as.formula(paste(rows, "~", cols)),
      scales = "free_x",
      labeller = ggplot2::label_both
    )
  p
}

#' Min-max score with 1 = best for a lower-is-better raw metric
#'
#' @noRd
.scale_lower_better <- function(x) {
  x <- as.numeric(x)
  finite <- x[is.finite(x)]
  if (length(finite) == 0L) {
    return(rep(NA_real_, length(x)))
  }
  lo <- min(finite)
  hi <- max(finite)
  rng <- hi - lo
  if (!is.finite(rng) || rng < .Machine$double.eps) {
    out <- rep(1, length(x))
    out[!is.finite(x)] <- NA_real_
    return(out)
  }
  1 - (x - lo) / rng
}

#' Pearson correlation of algorithms on Monte Carlo estimates
#'
#' @noRd
.one_algorithm_cor <- function(part, algos) {
  wide <- tidyr::pivot_wider(
    part,
    id_cols = dplyr::all_of(c("sample_id", "cell_type")),
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
  mat <- as.matrix(wide[, present, drop = FALSE])
  storage.mode(mat) <- "double"
  colnames(mat) <- present
  r_mat <- stats::cor(mat, use = "pairwise.complete.obs")
  .as_square_numeric_matrix(r_mat, labels = present)
}

#' Coerce a correlation table to a labelled square matrix
#'
#' @noRd
.as_square_numeric_matrix <- function(r_mat, labels = NULL) {
  r_mat <- as.matrix(unclass(r_mat))
  storage.mode(r_mat) <- "double"
  n <- nrow(r_mat)
  if (!identical(n, ncol(r_mat)) || n < 2L) {
    stop(
      "Correlation matrix must be square with at least two algorithms.",
      call. = FALSE
    )
  }
  if (is.null(labels)) {
    labels <- rownames(r_mat)
  }
  if (is.null(labels) || length(labels) != n) {
    labels <- colnames(r_mat)
  }
  if (is.null(labels) || length(labels) != n) {
    labels <- as.character(seq_len(n))
  }
  dimnames(r_mat) <- list(labels, labels)
  r_mat
}

#' Dissimilarity 1 - r for average-linkage clustering
#'
#' @noRd
.corr_distance <- function(r_mat) {
  r_mat <- .as_square_numeric_matrix(r_mat)
  d <- 1 - r_mat
  d[d < 0] <- 0
  diag(d) <- 0
  d
}

#' Long form of a named correlation matrix
#'
#' @noRd
.cor_to_tibble <- function(r_mat) {
  r_mat <- .as_square_numeric_matrix(r_mat)
  tbl <- as.data.frame(as.table(r_mat), stringsAsFactors = FALSE)
  names(tbl) <- c("algorithm_x", "algorithm_y", "correlation")
  tibble::as_tibble(tbl)
}

#' Hierarchical order from 1 - r
#'
#' @noRd
.hclust_corr_order <- function(r_mat) {
  d <- .corr_distance(r_mat)
  hc <- stats::hclust(stats::as.dist(d), method = "average")
  hc$labels[hc$order]
}

#' Algorithm-similarity correlation from a Monte Carlo benchmark
#'
#' Pearson correlations
#' \eqn{r_{ab}=\mathrm{cor}(\hat p_a,\hat p_b)} across Monte Carlo
#' replicates (and cell types). This is **behavioural similarity**: two
#' solvers can correlate near 1 while remaining systematically biased.
#' Hierarchical order uses \(d_{ab}=1-r_{ab}\).
#'
#' @param benchmark List from [run_simulation_benchmark()].
#' @param facet_rows,facet_cols Optional scenario columns. When supplied,
#'   correlations are computed inside each facet cell.
#'
#' @return A tibble with `algorithm_x`, `algorithm_y`, `correlation`,
#'   and any facet columns.
#'
#' @examples
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
  algos <- sort(unique(as.character(long$algorithm)))
  if (length(algos) < 2L) {
    stop(
      "algorithm_similarity() needs at least two algorithms.",
      call. = FALSE
    )
  }
  group_cols <- unique(c(facet_rows, facet_cols))
  group_cols <- group_cols[group_cols %in% names(long)]
  if (length(group_cols) == 0L) {
    r_mat <- .one_algorithm_cor(long, algos)
    return(.cor_to_tibble(r_mat))
  }
  keys <- interaction(long[, group_cols, drop = FALSE], drop = TRUE)
  pieces <- lapply(split(long, keys), function(part) {
    r_mat <- .one_algorithm_cor(part, algos)
    tbl <- .cor_to_tibble(r_mat)
    meta <- part[1L, group_cols, drop = FALSE]
    dplyr::bind_cols(meta, tbl)
  })
  dplyr::bind_rows(pieces)
}

#' Tile heatmap of algorithm-similarity correlations
#'
#' `ggplot2::geom_tile()` display of [algorithm_similarity()], with
#' rows and columns ordered by average-linkage clustering of
#' \(1-r\). Optional dendrogram via `ggdendro` (Suggests). This is the
#' default for a small correlation matrix; [plot_correlation_Heatmap()]
#' is reserved for linked multi-omics grids.
#'
#' @inheritParams algorithm_similarity
#' @param dendrogram If `TRUE`, attach a `ggdendro` ggplot as attribute
#'   `"dendrogram"` (ignored when scenario facets are used).
#'
#' @return A `ggplot` object.
#'
#' @examples
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
  .check_plot_dependencies(need_ggdist = FALSE)
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
  r_mat <- stats::xtabs(
    correlation ~ algorithm_x + algorithm_y,
    data = order_src
  )
  r_mat <- .as_square_numeric_matrix(r_mat)
  ord <- .hclust_corr_order(r_mat)
  sim$algorithm_x <- factor(sim$algorithm_x, levels = ord)
  sim$algorithm_y <- factor(sim$algorithm_y, levels = rev(ord))
  p <- ggplot2::ggplot(
    sim,
    ggplot2::aes(
      x = .data[["algorithm_x"]],
      y = .data[["algorithm_y"]],
      fill = .data[["correlation"]]
    )
  ) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::coord_equal() +
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
      fill = "Pearson r",
      caption = paste(
        "Clustering uses 1 - r (behavioural similarity),",
        "not numerical agreement of estimates."
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
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
      .check_plot_dependencies(need_ggdendro = TRUE)
      hc <- stats::hclust(
        stats::as.dist(.corr_distance(r_mat)),
        method = "average"
      )
      attr(p, "dendrogram") <- ggdendro::ggdendrogram(hc, rotate = TRUE)
    }
  }
  p
}

#' Long ADEMP table for faceted metric dots
#'
#' @noRd
.mc_metric_dot_table <- function(benchmark, metrics) {
  allowed <- c(
    "bias",
    "rmse",
    "mae",
    "coverage",
    "mean_interval_width"
  )
  metrics <- vapply(
    metrics,
    function(m) .match_arg_case_insensitive(m, allowed),
    character(1)
  )
  mc <- benchmark$monte_carlo
  long_est <- pivot_mc_estimates(benchmark)
  keys <- unique(c(
    "algorithm",
    "cell_type",
    .benchmark_meta_keys(benchmark$config, benchmark$optimisation)
  ))
  keys <- keys[keys %in% names(long_est)]
  mae_tbl <- long_est |>
    dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
    dplyr::summarise(
      mae = mean(abs(.data[["error"]]), na.rm = TRUE),
      .groups = "drop"
    )
  forest <- dplyr::left_join(mc, mae_tbl, by = intersect(names(mc), keys))
  pieces <- list()
  if ("bias" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "bias",
        raw = abs(.data[["bias"]])
      )
  }
  if ("rmse" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "rmse",
        raw = .data[["rmse"]]
      )
  }
  if ("mae" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "mae",
        raw = .data[["mae"]]
      )
  }
  if ("coverage" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "coverage",
        raw = abs(.data[["coverage"]] - 0.95)
      )
  }
  if ("mean_interval_width" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "mean_interval_width",
        raw = .data[["mean_interval_width"]]
      )
  }
  keep <- unique(c(keys, "metric", "raw"))
  dplyr::bind_rows(
    purrr::map(pieces, function(x) dplyr::select(x, dplyr::any_of(keep)))
  )
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
#' @examples
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
  .check_plot_dependencies(need_ggdist = FALSE)
  if (is.data.frame(benchmark)) {
    stop(
      "plot_mc_metric_dots() needs the full run_simulation_benchmark() list.",
      call. = FALSE
    )
  }
  plot_df <- .mc_metric_dot_table(benchmark, metrics)
  scale_cols <- unique(c("metric", facet_rows, facet_cols))
  scale_cols <- scale_cols[scale_cols %in% names(plot_df)]
  plot_df <- plot_df |>
    dplyr::group_by(dplyr::across(dplyr::all_of(scale_cols))) |>
    dplyr::mutate(score = .scale_lower_better(.data[["raw"]])) |>
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
    w_map <- weights
    comp <- plot_df |>
      dplyr::filter(.data[["metric"]] %in% w_names) |>
      dplyr::mutate(w = unname(w_map[as.character(.data[["metric"]])])) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(
        metric = "composite",
        raw = NA_real_,
        score = sum(.data[["w"]] * .data[["score"]], na.rm = TRUE),
        .groups = "drop"
      )
    plot_df <- dplyr::bind_rows(plot_df, comp)
  }
  plot_df$metric <- factor(
    plot_df$metric,
    levels = unique(c(metrics, if (!is.null(weights)) "composite"))
  )
  p <- ggplot2::ggplot(
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
    ggplot2::theme_bw()
  rows <- if (is.null(facet_rows)) "." else facet_rows
  cols <- if (is.null(facet_cols)) {
    "metric"
  } else {
    paste("metric", "+", facet_cols)
  }
  p +
    ggplot2::facet_grid(
      stats::as.formula(paste(rows, "~", cols)),
      labeller = ggplot2::label_both
    )
}
