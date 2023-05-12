test_that("test visualisations", {
  overlapping_unbalanced <- readRDS(test_path("fixtures", "temp_bivariate_1.rds"))
  overlapping_unbalanced <- overlapping_unbalanced %>% dplyr::select(-c(proportions, true_parameters, 
                                                                        overlap, variance, entropy))
  
  heatmap_overlapping_unbalanced <- plot_correlation_Heatmap (overlapping_unbalanced)

})
