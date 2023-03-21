# score_variable <- match.arg(score_variable, c("model_mse", "model_rmse", "model_coef_determination", "model_coef_determination_adjusted",
#                                               "model_mae", "model_cor", "model_ccc"))
# simulations_dp <- simulations_dp %>% dplyr::bind_rows(
#   tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
#                  variance="homoscedastic", overlap=overlap_homoscedastic, hellinger=hellinger_homoscedastic, density_plot=list(dp)))
#
# dp <-  plot_multivariate_density(simulated_data$X) +
#   ggtitle(paste(title_dp, "in hetero scenario"), subtitle = paste("overlap:", overlap_heteroscedastic,
#                                                                   "hellinger:", hellinger_heteroscedastic))
#
#
# simulations_dp <- simulations_dp %>% dplyr::bind_rows(
#   tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
#                  variance="heteroscedastic", overlap=overlap_heteroscedastic, hellinger=hellinger_heteroscedastic, density_plot=list(dp)))
#
# ##################################################################
# ##                       Complex Heatmaps                       ##
# ##################################################################
# CH_homo <- purrr::imap(split(simulation_metrics %>% filter(variance=="homoscedastic"), simulation_metrics$deconvolution_name), function(.x, .y) {
#   cor_matrix_per_algo <- .x %>% select (all_of(c("correlation_celltype1", "correlation_celltype2", score_variable))) %>%
#     tidyr::pivot_wider(names_from = c(correlation_celltype2), values_from = .data[[score_variable]],values_fn = mean) %>%
#     tibble::column_to_rownames("correlation_celltype1") %>% as.matrix()
#
#   complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(cor_matrix_per_algo, name=gsub("model_", "", score_variable),
#                                                       heatmap_legend_param = list(title = gsub("model_", "", score_variable)), row_title="Correlation cell type 1",
#                                                       cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix_per_algo),
#                                                       row_title_gp = grid::gpar(fontsize = 10), column_names_rot = 45,
#                                                       cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
#                                                       column_labels = colnames(cor_matrix_per_algo),
#                                                       width = unit(8, "cm"), height = unit(8, "cm"),
#                                                       column_title_gp = grid::gpar(fontsize = 10), column_title="Correlation cell type 2") %>%
#     ComplexHeatmap::draw(padding=unit(c(0, 0, 0, 0), "cm"),
#                          column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>%
#     grid::grid.grabExpr()
#   return(complex_heatmap_per_algo)
# })
# CH_homo <- gridExtra::arrangeGrob(grobs = CH_homo, ncol = min(length(deconvolution_functions), 2),  padding = unit(0.1, "line"),
#                                   top=ggpubr::text_grob(paste("Entropy is:", entropy, "\n homo scenario"), size = 18, face = "bold"))
#
# CH_hetero <- purrr::imap(split(simulation_metrics %>% filter(variance=="heteroscedastic"), simulation_metrics$deconvolution_name), function(.x, .y) {
#   cor_matrix_per_algo <- .x %>% select (all_of(c("correlation_celltype1", "correlation_celltype2", score_variable))) %>%
#     tidyr::pivot_wider(names_from = c(correlation_celltype2), values_from = .data[[score_variable]],values_fn = mean) %>%
#     tibble::column_to_rownames("correlation_celltype1") %>% as.matrix()
#
#   complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(cor_matrix_per_algo, name=gsub("model_", "", score_variable),
#                                                       heatmap_legend_param = list(title = gsub("model_", "", score_variable)), row_title="Correlation cell type 1",
#                                                       cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix_per_algo),
#                                                       row_title_gp = grid::gpar(fontsize = 10), column_names_rot = 45,
#                                                       cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
#                                                       column_labels = colnames(cor_matrix_per_algo),
#                                                       width = unit(8, "cm"), height = unit(8, "cm"),
#                                                       column_title_gp = grid::gpar(fontsize = 10), column_title="Correlation cell type 2") %>%
#     ComplexHeatmap::draw(padding=unit(c(0, 0, 0, 0), "cm"),
#                          column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>%
#     grid::grid.grabExpr()
#   return(complex_heatmap_per_algo)
# })
# CH_hetero <- gridExtra::arrangeGrob(grobs = CH_hetero, ncol = min(length(deconvolution_functions), 2),  padding = unit(0.1, "line"),
#                                     top=ggpubr::text_grob(paste("Entropy is:", entropy, "\n hetero scenario"), size = 18, face = "bold"))
#
#
# #################################################################
# ##                        mse boxplots                        ##
# #################################################################
#
# mse_boxplots <- plot_boxplots_parameters(simulation_metrics, score_variable) +
#   ggtitle(paste("Entropy is:", entropy))
#
# metric_plots <- gridExtra::arrangeGrob(grobs=list(CH_homo, CH_hetero, mse_boxplots),
#                                        layout_matrix = matrix(c(1, 2, 3, 3), nrow = 2, byrow = T), heights = c(1.5, 2))
