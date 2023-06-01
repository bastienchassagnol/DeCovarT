test_that("test visualisations", {
  bivariate_estimation <- readRDS(testthat::test_path("fixtures", "bivariate_estimation.rds"))
  heatmap_overlapping_unbalanced <- plot_correlation_Heatmap (bivariate_estimation %>% dplyr::filter(ID=="B1_Ho"))
  
  vdiffr::expect_doppelganger("Correlation heatmap plot", heatmap_overlapping_unbalanced$lm)
})
