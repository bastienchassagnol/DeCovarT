###########################################################################
###########################################################################
###                                                                     ###
###                   GENERATE HEATMAP VISUALISATIONS                   ###
###                                                                     ###
###########################################################################
###########################################################################

##################################################################
##      generate WABi general complexHeatmap       ##
##################################################################

library(dplyr)
bivariate_simulation <- readRDS("./simulations/results/bivariate_scenario.rds")
bivariate_simulation_test <- bivariate_simulation %>%
  mutate(algorithm=factor(algorithm, levels = c("lm", "nnls", "lsei",
                                                "gradient", "hessian", "DeCoVarT",
                                                "optim", "barrier", "SA"))) %>%
  mutate(algorithm=forcats::fct_recode(algorithm, Levenberg = "DeCoVarT") %>% 
           forcats::fct_relevel())
saveRDS(bivariate_simulation_test %>% dplyr::select(-c("proportions", "true_parameters",
                                                       "variance", "overlap", "model_coef_determination", 
                                                       "model_coef_determination_adjusted", "model_cor",
                                                       "entropy", "centroids")) %>% 
          dplyr::relocate("ID", .before = "correlation_celltype1") %>% 
          dplyr::rename_with(~gsub("celltype_", "p", .x)), 
        file = "./data/bivariate/bivariate_parameters.rds")

# bivariate_simulation <- bivariate_simulation %>%
#   mutate(ID = case_when(
#     proportions == "balanced" & centroids=="small CLD" ~ "B1_",
#     proportions == "highly_unbalanced" & centroids=="small CLD" ~ "B2_",
#     proportions == "balanced" & centroids=="high CLD" ~ "B3_",
#     proportions == "highly_unbalanced" & centroids=="high CLD" ~ "B4_",
#     TRUE                      ~ as.character(ID) )) %>% 
#   dplyr::mutate(ID = dplyr::if_else(variance=="homoscedasctic", paste0(ID, "Ho"), paste0(ID, "He")))
# saveRDS(bivariate_simulation, "./simulations/results/bivariate_scenario.rds")
# saveRDS(bivariate_formatted %>% dplyr::mutate(across(where(is.numeric), signif, 4)) %>% 
#           dplyr::select(-c("model_rmse", "model_mae")),
#         "./data/bivariate/bivariate_parameters.rds")

splitted_heatmap <- split(x = bivariate_simulation_test, f = bivariate_simulation_test$ID)
bivariate_simulation_heatmap <- purrr::imap(splitted_heatmap, function (.data, .name_scenario) {
  print(paste("Name scenario is ", .name_scenario))
  heatmap_per_scenario <- plot_correlation_Heatmap(.data)
  heatmap_page <- purrr::imap(heatmap_per_scenario,   ~ ComplexHeatmap::draw(.x, padding=unit(c(0, 0, 0, 0), "cm"),
                                                                     column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>%
                                grid::grid.grabExpr())
  heatmap_page <- gridExtra::arrangeGrob(grobs = heatmap_page, ncol = 3,  padding = unit(0.1, "line"),
                                         top=ggpubr::text_grob(.name_scenario, size = 18, face = "bold"))
  return(heatmap_page)
})

ggsave("./figs/bivariate_Heatmaps_test.pdf", 
       gridExtra::marrangeGrob(grobs=bivariate_simulation_heatmap, top="", ncol = 1, nrow = 1),
       width = 12, height = 12,dpi = 300)

  







  