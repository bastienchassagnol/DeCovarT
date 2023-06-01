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

test_bivariate_scenario <- withr::with_seed(seed = 3L,
                                            benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.50, 0.50)),
                                                                                         n = 2, corr_sequence = c(-0.75, 0.75),
                                                                                         signature_matrices=list("small CLD" = matrix(c(20, 22, 22, 20), nrow = 2)),
                                                                                         deconvolution_functions = deconvolution_functions),
                                            .rng_kind = "L'Ecuyer-CMRG")

saveRDS(test_bivariate_scenario$config, testthat::test_path("fixtures", "bivariate_configuraton.rds"))
saveRDS(test_bivariate_scenario$simulations, testthat::test_path("fixtures", "bivariate_estimation.rds"))