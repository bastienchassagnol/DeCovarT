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



# test <- highly_overlapping_simulation %>% 
#   filter(deconvolution_name=="DeCoVarT") %>% 
#   group_by(proportions, correlation_celltype1, correlation_celltype2) %>% 
#   summarise(num_simulations=n()) %>% 
#   filter(num_simulations!=500) %>% 
#   arrange(num_simulations)



non_overlapping_simulation %>% group_by(deconvolution_name) %>% 
  summarise(num_simulations=n())

  