##################################################################
##                   two-component simulation                   ##
##################################################################


test_that("small simulation testing", {
  simulation_two_genes <- simulate_bulk_mixture (signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2,
                                                                        dimnames = list(paste0("genes_", 1:2), paste0("cell_type_", 1:2))),
                                                           cov_tensor=array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2),
                                                                            dimnames = list(paste0("genes_", 1:2), paste0("genes_", 1:2), paste0("cell_type_", 1:2))), n=10)


})


test_that("understand mistakes", {
  erreur_output <- readRDS("./simulations/erreurs/erreur_1_function_DeCoVarT.rds")
  wrong_estimate <- do.call(deconvolute_ratios_DeCoVarT, erreur_output)
})


test_that("benchmark deconvolution", {
  RNGkind("L'Ecuyer-CMRG"); set.seed(3) # set an unique core
  deconvolution_functions <- list("gradient"=list(FUN=deconvolute_ratios_first_order, 
                                                  additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "lsei"=list(FUN=deconvolute_ratios_deconRNASeq),
                                  "lm"=list(FUN=deconvolute_ratios_abbas),
                                  "nnls"=list(FUN=deconvolute_ratios_nnls),
                                  # with the new log-likelihood function
                                  "optim"=list(FUN=deconvolute_ratios_basic_optim, 
                                               additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "barrier"=list(FUN=deconvolute_ratios_constrOptim,
                                                 additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  # "LBFGS"=list(FUN=deconvolute_ratios_LBFGS,
                                  #              additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  
                                  "gradient"=list(FUN=deconvolute_ratios_first_order, 
                                                  additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "hessian"=list(FUN=deconvolute_ratios_second_order,
                                                 additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT,
                                                  additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "SA"=list(FUN=deconvolute_ratios_simulated_annealing,
                                            additionnal_parameters = list(epsilon = 10^-3, itmax = 200)))
  
  easy_scenario <- benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.50, 0.50)), n = 1, 
                                                                corr_sequence = c(-0.8, -0.6),
                                                                signature_matrices=list("small CLD" = matrix(c(20, 22, 22, 20), nrow = 2),
                                                                                        "high CLD" = matrix(c(20, 40, 40, 20), nrow = 2)),
                                                                deconvolution_functions = deconvolution_functions)
  

})


test_that("benchmark metrics and plots", {
  # aggregate metrics
  mean_sd <- list(mean = ~mean(.x, na.rm = TRUE),
                  sd = ~sd(.x, na.rm = TRUE))
  simulation_metrics_summary <- two_genes_heteroscedastic_high_OVL %>% purrr::map_dfr("simulation_metrics")  %>%
    group_by(proportions, correlation_celltype1, correlation_celltype2, variance, deconvolution_name) %>%
    summarise(across(c(model_mse, model_rmse, model_coef_determination, model_mae, celltype_1, celltype_2), mean_sd)) %>%
    arrange(correlation_celltype1, correlation_celltype2, desc(proportions), desc(variance))
})



