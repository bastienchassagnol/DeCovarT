###########################################################################
###########################################################################
###                                                                     ###
###                   GENERATE HEATMAP VISUALISATIONS                   ###
###                                                                     ###
###########################################################################
###########################################################################


##################################################################
##          compare own deconvolution in easy scenario          ##
##################################################################


easy_simulation <- readRDS("./simulations/results/easy_scenario.rds")

easy_simulation_heatmap <- list()
for (scen_prop in unique(easy_simulation$proportions)) {
  for (scen_var in unique(easy_simulation$variance)) {
    name_scenario <- paste0(scen_var, " with ", scen_prop, " ratios") %>% stringr::str_to_title()
    heatmap_list <- easy_simulation %>% 
      dplyr::filter(proportions==scen_prop & variance==scen_var) %>% 
      plot_correlation_Heatmap()
    heatmap_page <- purrr::imap(heatmap_list,   ~ ComplexHeatmap::draw(.x, padding=unit(c(0, 0, 0, 0), "cm"),
                                                                       column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>% 
                                  grid::grid.grabExpr())
    easy_simulation_heatmap[[name_scenario]] <- gridExtra::arrangeGrob(grobs = heatmap_page, ncol = 2,  padding = unit(0.1, "line"),
                                                                      top=ggpubr::text_grob(name_scenario, size = 18, face = "bold"))
  }
}

ggsave("./figs/Heatmap_easy_example.pdf", 
       gridExtra::marrangeGrob(grobs=easy_simulation_heatmap, top="", ncol = 1, nrow = 1), width = 12, height = 20,dpi = 300)



#################################################################
##        compare own deconvolution in complex scenario        ##
#################################################################
complex_simulation <- readRDS("./simulations/results/complex_scenario.rds")

complex_simulation_heatmap <- list()
for (scen_prop in unique(complex_simulation$proportions)) {
  for (scen_var in unique(complex_simulation$variance)) {
    name_scenario <- paste0(scen_var, " with ", scen_prop, " ratios") %>% stringr::str_to_title()
    heatmap_list <- complex_simulation %>% 
      dplyr::filter(proportions==scen_prop & variance==scen_var) %>% 
      plot_correlation_Heatmap()
    heatmap_page <- purrr::imap(heatmap_list,   ~ ComplexHeatmap::draw(.x, padding=unit(c(0, 0, 0, 0), "cm"),
                                                                       column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>% 
                                  grid::grid.grabExpr())
    complex_simulation_heatmap[[name_scenario]] <- gridExtra::arrangeGrob(grobs = heatmap_page, ncol = 2,  padding = unit(0.1, "line"),
                                                                       top=ggpubr::text_grob(name_scenario, size = 18, face = "bold"))
  }
}

ggsave("./figs/Heatmap_complex_example.pdf", 
       gridExtra::marrangeGrob(grobs=complex_simulation_heatmap, top="", ncol = 1, nrow = 1), width = 12, height = 20,dpi = 300)


##################################################################
##       the true benchmarks, in non-overlapping scenario       ##
##################################################################


non_overlapping_simulation <- readRDS("./simulations/results/small_overlap_version2.rds")

non_overlapping_simulation_heatmap <- list()
for (scen_prop in unique(non_overlapping_simulation$proportions)) {
  for (scen_var in unique(non_overlapping_simulation$variance)) {
    name_scenario <- paste0(scen_var, " with ", scen_prop, " ratios") %>% stringr::str_to_title()
    heatmap_list <- non_overlapping_simulation %>% 
      dplyr::filter(proportions==scen_prop & variance==scen_var) %>% 
      plot_correlation_Heatmap(score_variable = "model_mse")
    heatmap_page <- purrr::imap(heatmap_list,   ~ ComplexHeatmap::draw(.x, padding=unit(c(0, 0, 0, 0), "cm"),
                                                                       column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>% 
                                  grid::grid.grabExpr())
    non_overlapping_simulation_heatmap[[name_scenario]] <- gridExtra::arrangeGrob(grobs = heatmap_page, ncol = 2,  padding = unit(0.1, "line"),
                                                                          top=ggpubr::text_grob(name_scenario, size = 18, face = "bold"))
  }
}

ggsave("./figs/Heatmap_non_overlapping_example_2.pdf", 
       gridExtra::marrangeGrob(grobs=non_overlapping_simulation_heatmap, top="", ncol = 1, nrow = 1), width = 12, height = 20,dpi = 300)


##################################################################
##       the true benchmarks, in highly-overlapping scenario       ##
##################################################################


highly_overlapping_simulation <- readRDS("./simulations/results/high_overlap.rds")

highly_overlapping_simulation_heatmap <- list()
for (scen_prop in unique(highly_overlapping_simulation$proportions)) {
  for (scen_var in unique(highly_overlapping_simulation$variance)) {
    name_scenario <- paste0(scen_var, " with ", scen_prop, " ratios") %>% stringr::str_to_title()
    heatmap_list <- highly_overlapping_simulation %>% 
      dplyr::filter(proportions==scen_prop & variance==scen_var) %>% 
      plot_correlation_Heatmap(score_variable = "model_mse")
    heatmap_page <- purrr::imap(heatmap_list,   ~ ComplexHeatmap::draw(.x, padding=unit(c(0, 0, 0, 0), "cm"),
                                                                       column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>% 
                                  grid::grid.grabExpr())
    highly_overlapping_simulation_heatmap[[name_scenario]] <- gridExtra::arrangeGrob(grobs = heatmap_page, ncol = 2,  padding = unit(0.1, "line"),
                                                                                  top=ggpubr::text_grob(name_scenario, size = 18, face = "bold"))
  }
}

ggsave("./figs/Heatmap_highly_overlapping.pdf", 
       gridExtra::marrangeGrob(grobs=highly_overlapping_simulation_heatmap, top="", ncol = 1, nrow = 2), width = 12, height = 20,dpi = 300)


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

p2 <- grid::grid.grabExpr(ComplexHeatmap::draw(global_heatmap))
cowplot_p2 <- cowplot::plot_grid(p2, nrow = 1, ncol=1)
)
global_heatmap <- global_heatmap %>%  ComplexHeatmap::draw() %>% grid::grid.grabExpr()







# test <- highly_overlapping_simulation %>% 
#   filter(deconvolution_name=="DeCoVarT") %>% 
#   group_by(proportions, correlation_celltype1, correlation_celltype2) %>% 
#   summarise(num_simulations=n()) %>% 
#   filter(num_simulations!=500) %>% 
#   arrange(num_simulations)




  