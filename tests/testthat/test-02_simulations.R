##################################################################
##                   two-component simulation                   ##
##################################################################


test_that("small simulation testing", {
  simulation_two_genes <- withr::with_seed(seed = 3L,
                                           simulate_bulk_mixture (signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2,
                                                                                          dimnames = list(paste0("genes_", 1:2), paste0("cell_type_", 1:2))),
                                                                  cov_tensor=array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2),
                                                                                   dimnames = list(paste0("genes_", 1:2), paste0("genes_", 1:2),
                                                                                                   paste0("cell_type_", 1:2))), n=10), 
                                           .rng_kind = "L'Ecuyer-CMRG")
  expect_equal(simulation_two_genes$Y[,1:2], 
               matrix(c(30.16652, 29.72223, 30.31822, 29.41423), nrow = 2,
                      dimnames = list(paste0("genes_", 1:2), paste0("sample_", 1:2))), tolerance = 10^-3)
  

  
})


test_that("benchmark deconvolution", {
  deconvolution_functions <- list("lm"=list(FUN=deconvolute_ratios_abbas),
                                  "nnls"=list(FUN=deconvolute_ratios_nnls),
                                  "lsei"=list(FUN=deconvolute_ratios_deconRNASeq),
                                  "LBFGS"=list(FUN=deconvolute_ratios_LBFGS,
                                               additional_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "gradient"=list(FUN=deconvolute_ratios_first_order,
                                                  additional_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "hessian"=list(FUN=deconvolute_ratios_second_order,
                                                 additional_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT,
                                                  additional_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  # with the new log-likelihood function
                                  "optim"=list(FUN=deconvolute_ratios_basic_optim,
                                               additional_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "barrier"=list(FUN=deconvolute_ratios_constrOptim,
                                                 additional_parameters = list(epsilon = 10^-3, itmax = 200)),
                                  "SA"=list(FUN=deconvolute_ratios_simulated_annealing,
                                            additional_parameters = list(epsilon = 10^-3, itmax = 200)))
  
  bivariate_scenario <- withr::with_seed(seed = 3L,
                                         suppressMessages(benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.50, 0.50)),
                                                                                      n = 2, corr_sequence = c(-0.75, 0.75),
                                                                                      signature_matrices=list("small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)),
                                                                                      deconvolution_functions = deconvolution_functions)),
                                         .rng_kind = "L'Ecuyer-CMRG")
  bivariate_configuraton <- readRDS(testthat::test_path("fixtures", "bivariate_configuraton.rds"))
  bivariate_estimation <- readRDS(testthat::test_path("fixtures", "bivariate_estimation.rds"))
  
  expect_equal(bivariate_configuraton, bivariate_scenario$config)
  expect_equal(bivariate_estimation, bivariate_scenario$simulations)

})




