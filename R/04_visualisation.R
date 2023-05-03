#' Plot Correlation Heatmap
#'
#' For each deconvolution algorithm, plot the metric selected over the range of selected ratios
#'
#' @param distribution_metrics A tibble with the metric scores
#' @param score_variable The name of the metric to plot on the Heatmap
#'
#' @export

plot_correlation_Heatmap <- function(distribution_metrics, score_variable = "model_rmse") {
  score_variable <- match.arg(score_variable, c("model_mse", "model_rmse", "model_coef_determination",
                                                "model_coef_determination_adjusted",
                                                "model_mae", "model_cor", "model_ccc"))
  distribution_metrics <- distribution_metrics %>% 
    dplyr::select(dplyr::all_of(c("correlation_celltype1", "correlation_celltype2", "deconvolution_name", score_variable)))
  
  
  complex_heatmap_list <- purrr::imap(split(distribution_metrics, distribution_metrics[["deconvolution_name"]]), function(.x, .y) {
    cor_matrix_per_algo <- .x %>% dplyr::select(-"deconvolution_name") %>% 
      tidyr::pivot_wider(names_from = c(correlation_celltype2), 
                         values_from = .data[[score_variable]],values_fn = mean) %>%
      tibble::column_to_rownames("correlation_celltype1") %>% as.matrix()
    
    complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(cor_matrix_per_algo, name=gsub("model_", "", score_variable),
                                                        heatmap_legend_param = list(title = gsub("model_", "", score_variable) %>% toupper),
                                                        row_title="Corr cell type 1", cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8),
                                                        row_labels = colnames(cor_matrix_per_algo), row_title_gp = grid::gpar(fontsize = 10), 
                                                        column_names_rot = 0, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
                                                        column_labels = colnames(cor_matrix_per_algo),
                                                        width = unit(8, "cm"), height = unit(8, "cm"),
                                                        column_title_gp = grid::gpar(fontsize = 10), column_title="Corr cell type 2") 
    return(complex_heatmap_per_algo)
  })
  return(complex_heatmap_list)
}





