test_that("test visualisations", {
 
  toy_simulation <- readRDS("./simulations/results/toy_example_easy_scenario.rds")
  
  toy_simulation_heatmap <- list()
  for (scen_prop in unique(toy_simulation$proportions)) {
    for (scen_var in unique(toy_simulation$variance)) {
      name_scenario <- paste0(scen_var, " with ", scen_prop, " ratios") %>% stringr::str_to_title()
      heatmap_list <- toy_simulation %>% 
        dplyr::filter(proportions==scen_prop & variance==scen_var) %>% 
        plot_correlation_Heatmap()
      heatmap_page <- purrr::imap(heatmap_list,   ~ ComplexHeatmap::draw(.x, padding=unit(c(0, 0, 0, 0), "cm"),
                                                         column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>% 
                                    grid::grid.grabExpr())
      toy_simulation_heatmap[[name_scenario]] <- gridExtra::arrangeGrob(grobs = heatmap_page, ncol = 2,  padding = unit(0.1, "line"),
                                                                     top=ggpubr::text_grob(name_scenario, size = 18, face = "bold"))
    }
  }
  ggsave("./figs/Heatmap_toy_example.pdf", 
         gridExtra::arrangeGrob(grobs=toy_simulation_heatmap, top="", ncol = 1), width = 12, height = 15,dpi = 300)
  
  
})
