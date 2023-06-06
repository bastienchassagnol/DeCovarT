#' Plot Correlation Heatmap
#'
#' For each deconvolution algorithm, plot the metric selected over the range of selected ratios
#'
#' @param distribution_metrics A tibble with the metric scores
#' @param score_variable The name of the metric to be represented on the Heatmap
#' @param n_break the continous number of breaks allowed to generate the Heatmap
#' @param uni_scale if FALSE, each Heatmap is plotted with its own scale
#'
#' @export

plot_correlation_Heatmap <- function(distribution_metrics, score_variable = "model_mse",
                                     n_break = 20, uni_scale = TRUE) {
  score_variable <- match.arg(score_variable, c(
    "model_mse", "model_rmse", "model_coef_determination",
    "model_coef_determination_adjusted",
    "model_mae", "model_cor", "model_ccc"
  ))
  distribution_metrics <- distribution_metrics %>%
    dplyr::select(dplyr::all_of(c("correlation_celltype1", "correlation_celltype2", "algorithm", score_variable))) %>%
    dplyr::mutate(algorithm = factor(.data$algorithm, levels = unique(.data$algorithm)))

  # design the scaling colour
  mean_distribution_metrics <- distribution_metrics %>%
    dplyr::group_by(.data$correlation_celltype1, .data$correlation_celltype2, .data$algorithm) %>%
    dplyr::summarise(mean_metric = mean(.data[[score_variable]], na.rm = TRUE)) # we remove missing cases

  min_metric <- min(mean_distribution_metrics %>% dplyr::pull("mean_metric"), na.rm = TRUE)
  max_metric <- max(mean_distribution_metrics %>% dplyr::pull("mean_metric"), na.rm = TRUE)
  if (uni_scale) {
    col <- circlize::colorRamp2(
      seq(min_metric, max_metric, length.out = n_break),
      viridis::viridis(n_break)
    )
  }

  complex_heatmap_list <- purrr::imap(split(distribution_metrics, distribution_metrics[["algorithm"]]), function(.x, .y) {
    cor_matrix_per_algo <- .x %>%
      dplyr::select(-"algorithm") %>%
      tidyr::pivot_wider(
        names_from = dplyr::all_of("correlation_celltype2"),
        values_from = dplyr::all_of(score_variable), values_fn = mean
      ) %>%
      tibble::column_to_rownames("correlation_celltype1") %>%
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

    complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(cor_matrix_per_algo,
      col = col, name = gsub("model_", "", score_variable),
      heatmap_legend_param = list(title = gsub("model_", "", score_variable) %>% toupper()),
      row_title = "Corr cell type 1", cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8),
      row_labels = colnames(cor_matrix_per_algo), row_title_gp = grid::gpar(fontsize = 10),
      column_names_rot = 0, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
      column_title_gp = grid::gpar(fontsize = 10), column_title = "Corr cell type 2",
      column_labels = colnames(cor_matrix_per_algo),
      width = unit(6, "cm"), height = unit(6, "cm")
    )
    return(complex_heatmap_per_algo)
  })
  return(complex_heatmap_list)
}
