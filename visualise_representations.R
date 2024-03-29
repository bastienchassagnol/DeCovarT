###########################################################################
###########################################################################
###                                                                     ###
###                   GENERATE HEATMAP VISUALISATIONS                   ###
###                                                                     ###
###########################################################################
###########################################################################

#################################################################
##            generate parameter configuration file            ##
#################################################################

library(dplyr); library(kableExtra)

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

bivariate_simulation <- readRDS("./simulations/results/bivariate_scenario.rds")

bivariate_configuration <- bivariate_simulation %>% select(c("ID", "overlap", "entropy",
                                                             "proportions", "variance",  "centroids", 
                                                             "true_parameters")) %>% 
  dplyr::distinct() %>% dplyr::mutate(centroids=stringr::str_replace_all(centroids, "CLD", "ICD"))
saveRDS(bivariate_configuration, "./simulations/results/complete_bivariate_configuration.rds")

reduced_bivariate_configuration <- bivariate_configuration %>% 
  mutate(ID=factor(ID, levels= unique(bivariate_configuration$ID))) %>% group_by(ID) %>% 
  summarise(Entropy = entropy, OVL = mean(overlap),
    Proportions = purrr::map_chr(true_parameters, ~ paste(.x$p, collapse = " / ")),
    Means = purrr::map_chr(true_parameters, ~ paste0("(", paste0(.x$mu[, 1], collapse = ","), ");(", paste0(.x$mu[, 2], collapse = ","), ")")),
    Variance = purrr::map_chr(true_parameters, ~ paste(c(.x$sigma[1, 1, 1], .x$sigma[2, 2, 1]), collapse = " / "))) %>% 
  dplyr::distinct()
saveRDS(reduced_bivariate_configuration, "./simulations/results/reduced_bivariate_configuration.rds")

reduced_bivariate_configuration %>%  kbl(
    booktabs = T, caption = "The 8 general scenarios tested to compare the performance of DeCovarT
    vs standard linear deconvolution model", escape = F, align = "c") %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>%
  row_spec(0, bold = T) %>%
  row_spec(1:8, hline_after = T) %>%
  kable_styling(bootstrap_options = c("hover", "condensed"))

##################################################################
##          reduce size of the parameter boostrap file          ##
##################################################################

reduced_bivariate_simulation <- bivariate_simulation %>%
  mutate(algorithm=factor(algorithm, levels = c("lm", "nnls", "lsei",
                                                "gradient", "hessian", "DeCoVarT",
                                                "optim", "barrier", "SA"))) %>%
  mutate(algorithm=forcats::fct_recode(algorithm, Levenberg = "DeCoVarT") %>% 
           forcats::fct_relevel()) %>% 
  mutate(across(where(is.numeric), signif(digits = 4))) %>% 
  select(-c("proportions", "true_parameters",
            "variance", "overlap", "model_coef_determination", 
            "model_coef_determination_adjusted", "model_cor",
            "entropy", "centroids")) %>% 
  dplyr::relocate("ID", .before = "correlation_celltype1") %>% 
  dplyr::rename_with(~gsub("celltype_", "p", .x))

saveRDS(reduced_bivariate_simulation, file = "./data/bivariate/bivariate_parameters.rds")

##################################################################
##      generate WABi general complexHeatmap       ##
##################################################################
load("./data/bivariate_parameters.rda")

splitted_parameters <- split(x = bivariate_parameters, f = bivariate_parameters$ID)
bivariate_simulation_heatmap <- purrr::imap(splitted_parameters, function(.data, .name_scenario) {
  heatmap_per_scenario <- plot_correlation_Heatmap(.data) # actual call to the associated DeCovarT function
  heatmap_page <- purrr::imap(heatmap_per_scenario, ~ ComplexHeatmap::draw(.x,
                                                                           padding = unit(c(0, 0, 0, 0), "cm"),
                                                                           column_title = .y, column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
  ) %>%
    grid::grid.grabExpr())
  # general organisation: 3 deconvolution algorithms per column
  heatmap_page <- gridExtra::arrangeGrob(
    grobs = heatmap_page, ncol = 3, padding = unit(0.1, "line"),
    top = ggpubr::text_grob(.name_scenario, size = 18, face = "bold")
  )
  return(heatmap_page)
})

# save the actual output
ggsave("./figs/bivariate_Heatmaps_test.pdf",
       gridExtra::marrangeGrob(grobs = bivariate_simulation_heatmap, top = "", ncol = 1, nrow = 1),
       width = 12, height = 12, dpi = 300
)

#################################################################
##                    generate JOBIM figure                    ##
#################################################################

highly_overlapping_simulation <- readRDS("./simulations/results/high_overlap_version2.rds")


data <- highly_overlapping_simulation %>% 
  dplyr::filter(proportions=="balanced" & variance=="homoscedasctic") %>% 
  plot_correlation_Heatmap(score_variable = "model_mse")
data1 <- data$lsei@matrix; data2 <- data$`tricky optim`@matrix


common_min = min(c(data1, data2)); common_max = max(c(data1, data2))
col_fun = circlize::colorRamp2(c(common_min, common_max), c("blue", "red"))

global_heatmap <- ComplexHeatmap::Heatmap(data1, col=col_fun, 
                                          heatmap_legend_param = list(title = "MSE"),
                                          row_title="Corr cell type 1", cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8),
                                          row_labels = colnames(data1), row_title_gp = grid::gpar(fontsize = 10), 
                                          column_names_rot = 0, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
                                          column_labels = colnames(data1), width = unit(8, "cm"), height = unit(8, "cm"),
                                          column_title_gp = grid::gpar(fontsize = 10), column_title="Corr cell type 2") + 
  ComplexHeatmap::Heatmap(data2, col=col_fun, show_heatmap_legend = F,
                          heatmap_legend_param = list(title = "MSE"),
                          row_title="Corr cell type 1", cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8),
                          row_labels = colnames(data2), row_title_gp = grid::gpar(fontsize = 10), 
                          column_names_rot = 0, cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
                          column_labels = colnames(data2), width = unit(8, "cm"), height = unit(8, "cm"),
                          column_title_gp = grid::gpar(fontsize = 10), column_title="Corr cell type 2") 




  