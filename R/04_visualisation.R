#' Plot deconvolution metric heatmaps
#'
#' @description
#' For each algorithm, visualises a selected score over the design grid of
#' simulated \eqn{\boldsymbol{p}} (and related scenario factors).
#'
#' @param distribution_metrics Tibble of metric scores from a benchmark.
#' @param score_variable Column name of the metric to display.
#' @param n_break Number of colour breaks.
#' @param uni_scale If `FALSE`, each panel uses its own colour scale.
#'
#' @importFrom rlang .data
#' @export
plot_correlation_Heatmap <- function(
  distribution_metrics,
  score_variable = "model_mse",
  n_break = 20,
  uni_scale = TRUE
) {
  score_variable <- match.arg(
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

  min_metric <- min(
    mean_distribution_metrics |> dplyr::pull("mean_metric"),
    na.rm = TRUE
  )
  max_metric <- max(
    mean_distribution_metrics |> dplyr::pull("mean_metric"),
    na.rm = TRUE
  )
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

      complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(
        cor_matrix_per_algo,
        col = col,
        name = gsub("model_", "", score_variable),
        heatmap_legend_param = list(
          title = gsub("model_", "", score_variable) |> toupper()
        ),
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
  return(complex_heatmap_list)
}
