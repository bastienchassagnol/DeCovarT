# source(Rscript_HPC_API)
library(DeCovarT); library(dplyr)

# Rscript_HPC ("/home/bncl_cb/rstudio/working/DeCovarT/launch_simulations.R", JobName = "simulation_bivariate")

RNGkind("L'Ecuyer-CMRG"); set.seed(3) # set an unique core

deconvolution_functions <- list("lm"=list(FUN=deconvolute_ratios_abbas),
                                "nnls"=list(FUN=deconvolute_ratios_nnls),
                                "lsei"=list(FUN=deconvolute_ratios_deconRNASeq),
                                # with the new log-likelihood function
                                "optim"=list(FUN=deconvolute_ratios_basic_optim, 
                                             additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                "barrier"=list(FUN=deconvolute_ratios_constrOptim,
                                               additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                "SA"=list(FUN=deconvolute_ratios_simulated_annealing,
                                          additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                # "LBFGS"=list(FUN=deconvolute_ratios_LBFGS,
                                #              additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                
                                "gradient"=list(FUN=deconvolute_ratios_first_order, 
                                                additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                "hessian"=list(FUN=deconvolute_ratios_second_order,
                                               additionnal_parameters = list(epsilon = 10^-3, itmax = 200)),
                                "DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT,
                                                additionnal_parameters = list(epsilon = 10^-3, itmax = 200)))

test_bivariate_scenario <- benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.50, 0.50)),
                                                                        n = 1, corr_sequence = c(-0.75, 0.75),
                                                                        signature_matrices=list("small CLD" = matrix(c(20, 22, 22, 20), nrow = 2),
                                                                                                "high CLD" = matrix(c(20, 40, 40, 20), nrow = 2)),
                                                                        deconvolution_functions = deconvolution_functions[1:3])

bivariate_scenario <- suppressWarnings(benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.50, 0.50), 
                                                                                                     "highly_unbalanced"=c(0.95, 0.05)), n = 2000, 
                                                              corr_sequence = seq(-0.75, 0.75, 0.25),
                                                              signature_matrices=list("small CLD" = matrix(c(20, 22, 22, 20), nrow = 2),
                                                                                      "high CLD" = matrix(c(20, 40, 40, 20), nrow = 2)),
                                                              deconvolution_functions = deconvolution_functions))
saveRDS(bivariate_scenario, "/home/bncl_cb/rstudio/working/DeCovarT/simulations/results/bivariate_scenario.rds")










