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
    .ui_abort(c(
      "{.fn plot_correlation_Heatmap} requires {.pkg {missing}}.",
      "i" = "Install CRAN dependencies with install.packages().",
      "i" = "Install ComplexHeatmap with BiocManager::install(\"ComplexHeatmap\")."
    ))
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

#' Check ggplot2 (Imports) and optional ggdist / ggdendro Suggests
#'
#' @keywords internal
#' @noRd
.check_plot_dependencies <- function(
  need_ggdist = FALSE,
  need_ggdendro = FALSE,
  need_cowplot = FALSE
) {
  pkgs <- "ggplot2"
  if (isTRUE(need_ggdist)) {
    pkgs <- c(pkgs, "ggdist")
  }
  if (isTRUE(need_ggdendro)) {
    pkgs <- c(pkgs, "ggdendro")
  }
  if (isTRUE(need_cowplot)) {
    pkgs <- c(pkgs, "cowplot")
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

#' Canonical factor orders for fig02 scenario columns and solvers
#'
#' @keywords internal
#' @noRd
.scenario_level_orders <- function() {
  list(
    centroids = c("small_CLD", "large_CLD"),
    variance = c("homoscedastic", "heteroscedastic"),
    proportions = c(
      "balanced",
      "moderately unbalanced",
      "highly unbalanced"
    )
  )
}

#' Solver display order (missing levels are dropped)
#'
#' @keywords internal
#' @noRd
.algorithm_level_order <- function() {
  c(
    "nnls",
    "lsei",
    "SA",
    "gradient",
    "LBFGS",
    "LBFGSB",
    "L-BFGS-B",
    "Newton-Raphson",
    "Marquardt-Levenberg"
  )
}

#' Relevel a vector, keeping only levels that appear
#'
#' @keywords internal
#' @noRd
.relevel_existing <- function(x, order) {
  x_chr <- as.character(x)
  present <- intersect(order, unique(x_chr))
  extra <- setdiff(unique(x_chr), present)
  levels <- c(present, extra)
  if (length(levels) == 0L) {
    return(factor(x_chr))
  }
  if (requireNamespace("forcats", quietly = TRUE) && length(present) > 0L) {
    forcats::fct_relevel(factor(x_chr, levels = levels), present)
  } else {
    factor(x_chr, levels = levels)
  }
}

#' @keywords internal
#' @noRd
.relevel_algorithm <- function(x) {
  .relevel_existing(x, .algorithm_level_order())
}

#' Relevel scenario and algorithm columns used in fig02 plots
#'
#' @keywords internal
#' @noRd
.relevel_scenario_table <- function(tbl) {
  if (is.null(tbl) || !is.data.frame(tbl)) {
    return(tbl)
  }
  orders <- .scenario_level_orders()
  for (nm in names(orders)) {
    if (nm %in% names(tbl)) {
      tbl[[nm]] <- .relevel_existing(tbl[[nm]], orders[[nm]])
    }
  }
  if ("algorithm" %in% names(tbl)) {
    tbl$algorithm <- .relevel_algorithm(tbl$algorithm)
  }
  tbl
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
  .relevel_scenario_table(long)
}

#' Fill missing `theta_true` from a theta artefact table
#'
#' @noRd
.fill_theta_true_from_theta <- function(benchmark) {
  cfg <- benchmark$config
  truths <- benchmark$theta_true
  if (!is.null(truths) && length(truths) == nrow(cfg)) {
    return(benchmark)
  }
  th_tbl <- benchmark$theta
  if (is.null(th_tbl) || !"true_theta" %in% names(th_tbl)) {
    return(benchmark)
  }
  if ("ID" %in% names(cfg) && "ID" %in% names(th_tbl)) {
    idx <- match(as.character(cfg$ID), as.character(th_tbl$ID))
    benchmark$theta_true <- lapply(
      th_tbl$true_theta[idx],
      .unwrap_true_theta
    )
    return(benchmark)
  }
  if (nrow(th_tbl) == nrow(cfg)) {
    benchmark$theta_true <- lapply(th_tbl$true_theta, .unwrap_true_theta)
  }
  benchmark
}

#' Expected-Fisher Wald SE at the true composition of one scenario
#'
#' @noRd
.theoretical_se_from_theta <- function(true_theta) {
  th <- .unwrap_true_theta(true_theta)
  nms <- .truth_cell_names(th)
  p <- stats::setNames(as.numeric(th$p), nms)
  se <- .ilr_wald_se(p, th$mu, th$sigma, warn = FALSE)
  tibble::tibble(
    cell_type = nms,
    theoretical_se = as.numeric(se[nms])
  )
}

#' Refresh Wald coverage from expected Fisher at \eqn{p^{\star}}
#'
#' [expected_fisher_unconstrained()] / [vcov_ilr_delta()] depend only on
#' \eqn{(p,\mu,\Sigma)}, so the Cramer--Rao SE is constant for a
#' scenario (the same object as [confint.decovart_fit()]). Coverage is
#' the Monte Carlo rate at which
#' \eqn{\hat p_j\pm z\,\mathrm{SE}_j(p^{\star})} covers
#' \eqn{p_j^{\star}}. Mean-only solvers are left unchanged.
#'
#' @noRd
.attach_expected_fisher_wald <- function(benchmark, level = 0.95) {
  z <- stats::qnorm((1 + level) / 2)
  cfg <- benchmark$config
  mc <- benchmark$monte_carlo
  if (is.null(mc) || is.null(cfg) || nrow(mc) == 0L) {
    return(benchmark)
  }
  benchmark <- .fill_theta_true_from_theta(benchmark)
  truths <- benchmark$theta_true
  if (is.null(truths) || length(truths) == 0L) {
    return(benchmark)
  }
  se_pieces <- lapply(seq_along(truths), function(i) {
    row <- .theoretical_se_from_theta(truths[[i]])
    if ("ID" %in% names(cfg)) {
      row$ID <- as.character(cfg$ID[[i]])
    }
    row
  })
  se_tbl <- dplyr::bind_rows(se_pieces)
  long <- pivot_mc_estimates(benchmark)
  join_se <- intersect(c("ID", "cell_type"), names(long))
  join_se <- join_se[join_se %in% names(se_tbl)]
  if ("ID" %in% names(long)) {
    long$ID <- as.character(long$ID)
  }
  if ("ID" %in% names(se_tbl)) {
    se_tbl$ID <- as.character(se_tbl$ID)
  }
  long <- dplyr::left_join(long, se_tbl, by = join_se)
  long$is_wald <- .uses_ilr_wald_algorithm(long$algorithm)
  long$covered <- long$is_wald &
    is.finite(long$estimate) &
    is.finite(long$p_true) &
    is.finite(long$theoretical_se) &
    abs(long$estimate - long$p_true) <= z * long$theoretical_se
  method <- "wilson"
  if ("coverage_interval" %in% names(mc)) {
    methods <- unique(stats::na.omit(as.character(mc$coverage_interval)))
    if (length(methods) == 1L) {
      method <- methods[[1L]]
    }
  }
  grp <- intersect(c("ID", "algorithm", "cell_type"), names(long))
  wald_long <- long[long$is_wald, , drop = FALSE]
  if (nrow(wald_long) == 0L) {
    return(benchmark)
  }
  grouped <- dplyr::group_by(wald_long, dplyr::across(dplyr::all_of(grp)))
  upd <- dplyr::summarise(
    grouped,
    theoretical_se = mean(.data[["theoretical_se"]], na.rm = TRUE),
    covered_list = list(.data[["covered"]]),
    .groups = "drop"
  )
  intervals <- lapply(
    upd$covered_list,
    function(cv) coverage_mc_interval(cv, method = method)
  )
  upd$coverage <- vapply(intervals, `[[`, numeric(1), "coverage")
  upd$coverage_lower <- vapply(intervals, `[[`, numeric(1), "lower")
  upd$coverage_upper <- vapply(intervals, `[[`, numeric(1), "upper")
  upd$mcse_coverage <- vapply(intervals, `[[`, numeric(1), "mcse")
  upd$coverage_interval <- vapply(intervals, `[[`, character(1), "method")
  upd$mean_model_se <- upd$theoretical_se
  upd$mean_model_sd <- upd$theoretical_se
  upd$mean_interval_width <- 2 * z * upd$theoretical_se
  upd$covered_list <- NULL
  replace_cols <- c(
    "theoretical_se",
    "coverage",
    "coverage_lower",
    "coverage_upper",
    "coverage_interval",
    "mcse_coverage",
    "mean_model_se",
    "mean_model_sd",
    "mean_interval_width"
  )
  is_wald_mc <- .uses_ilr_wald_algorithm(mc$algorithm)
  mc_wald <- mc[is_wald_mc, , drop = FALSE]
  mc_rest <- mc[!is_wald_mc, , drop = FALSE]
  drop_now <- intersect(replace_cols, names(mc_wald))
  if (length(drop_now) > 0L) {
    mc_wald <- mc_wald[, setdiff(names(mc_wald), drop_now), drop = FALSE]
  }
  join_mc <- intersect(grp, names(mc_wald))
  if ("ID" %in% names(mc_wald)) {
    mc_wald$ID <- as.character(mc_wald$ID)
  }
  if ("ID" %in% names(upd)) {
    upd$ID <- as.character(upd$ID)
  }
  mc_wald <- dplyr::left_join(mc_wald, upd, by = join_mc)
  if ("empirical_sd" %in% names(mc_wald)) {
    mc_wald$se_sd_ratio <- mc_wald$theoretical_se / mc_wald$empirical_sd
  }
  if (!"theoretical_se" %in% names(mc_rest)) {
    mc_rest$theoretical_se <- NA_real_
  }
  benchmark$monte_carlo <- dplyr::bind_rows(mc_wald, mc_rest)
  benchmark$theta_true <- truths
  benchmark
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
#' @seealso [plot_mc_forest()], [plot_mc_raincloud()],
#'   [plot_bivariate_metric_tiles()]
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

#' @noRd
.collapse_facet_vars <- function(x) {
  if (is.null(x) || length(x) == 0L) {
    return(".")
  }
  paste(x, collapse = " + ")
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
  for (nm in vars) {
    if (!nm %in% names(data)) {
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
  ggplot2::facet_grid(
    stats::as.formula(
      paste(
        .collapse_facet_vars(facet_rows),
        "~",
        .collapse_facet_vars(facet_cols)
      )
    ),
    labeller = ggplot2::label_value,
    scales = scales,
    drop = FALSE
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
  slab_alpha = 1
) {
  .check_plot_dependencies(need_ggdist = TRUE)
  quantity <- .match_arg_case_insensitive(quantity, c("error", "estimate"))
  df <- if (is.data.frame(benchmark)) {
    .relevel_scenario_table(benchmark)
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
  x_var <- if (identical(quantity, "error")) "error" else "estimate"
  x_lab <- if (identical(quantity, "error")) {
    "Monte Carlo error (estimate minus truth)"
  } else {
    "Monte Carlo estimate"
  }
  dodge <- ggplot2::position_dodge(width = dodge_width)
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
      justification = -0.12,
      point_interval = ggdist::median_qi,
      normalize = "groups",
      scale = slab_scale,
      interval_size = 2.8,
      point_size = 1.6,
      alpha = slab_alpha,
      position = dodge
    )
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
  if (identical(quantity, "error")) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
  } else {
    truth <- dplyr::distinct(
      df,
      .data[["cell_type"]],
      .data[["p_true"]]
    )
    ct_lvls <- unique(as.character(truth$cell_type))
    pal <- c("#E41A1C", "#4DAF4A", "#377EB8", "#984EA3")
    pal <- pal[seq_len(length(ct_lvls))]
    names(pal) <- ct_lvls
    for (ct in ct_lvls) {
      xs <- unique(truth$p_true[as.character(truth$cell_type) == ct])
      p <- p +
        ggplot2::geom_vline(
          xintercept = xs,
          colour = unname(pal[[ct]]),
          linetype = "longdash",
          linewidth = 0.7
        )
    }
    p <- p +
      ggplot2::labs(
        caption = paste(
          "Central 50% and 95% of Monte Carlo replicates;",
          "not a confidence interval for p.",
          "Dashed vertical lines: true cell-type proportions",
          "(type 1 red, type 2 green)."
        )
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
  .check_plot_dependencies(need_ggdist = FALSE)
  allowed <- c(
    "bias",
    "rmse",
    "mae",
    "coverage",
    "mean_interval_width",
    "se_sd_ratio",
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
  if ("se_sd_ratio" %in% metrics) {
    pieces[[length(pieces) + 1L]] <- forest |>
      dplyr::mutate(
        metric = "se_sd_ratio",
        estimate = .data[["se_sd_ratio"]],
        lower = NA_real_,
        upper = NA_real_,
        reference = 1
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
  plot_df <- dplyr::filter(plot_df, is.finite(.data[["estimate"]]))
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
  id_cols <- c("sample_id", "cell_type")
  if ("ID" %in% names(part)) {
    id_cols <- c("ID", id_cols)
  }
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
#' Hierarchical order uses \eqn{d_{ab} = 1 - r_{ab}}.
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
  long$algorithm <- .relevel_algorithm(long$algorithm)
  algos <- levels(droplevels(long$algorithm))
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
#' \(1-r\). Optional dendrogram via `ggdendro` (Suggests), drawn to
#' the **right** of the tiles with leaves flush against the heatmap.
#' This is the default for a small correlation matrix;
#' [plot_correlation_Heatmap()] is reserved for linked multi-omics grids.
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
      fill = "Pearson r",
      caption = paste(
        "Clustering uses 1 - r (behavioural similarity),",
        "not numerical agreement of estimates."
      )
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
      .check_plot_dependencies(
        need_ggdendro = TRUE,
        need_cowplot = TRUE
      )
      hc <- stats::hclust(
        stats::as.dist(.corr_distance(r_mat)),
        method = "average"
      )
      dend <- .similarity_dendrogram_plot(hc, n_leaf = length(ord))
      tiles <- p +
        ggplot2::theme(
          plot.margin = ggplot2::margin(4, 0, 4, 4)
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
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      limits = c(0.5, n_leaf + 0.5),
      expand = c(0, 0)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        colour = NA
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
      plot.margin = ggplot2::margin(4, 4, 4, 0)
    )
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
    theme_decovart_facets()
  p +
    .facet_grid_from_names(
      plot_df,
      facet_rows,
      unique(c("metric", facet_cols))
    )
}
