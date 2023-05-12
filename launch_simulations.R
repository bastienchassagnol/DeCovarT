# source(Rscript_HPC_API)
library(DeCovarT); library(dplyr)

# Rscript_HPC ("/home/bncl_cb/rstudio/working/DeCovarT/launch_simulations.R", JobName = "simulation_bivariate")

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

bivariate_scenario <- suppressWarnings(benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.50, 0.50), "highly_unbalanced"=c(0.95, 0.05)), n = 2000, 
                                                              corr_sequence = seq(-0.75, 0.75, 0.25),
                                                              signature_matrices=list("small CLD" = matrix(c(20, 22, 22, 20), nrow = 2),
                                                                                      "high CLD" = matrix(c(20, 40, 40, 20), nrow = 2)),
                                                              deconvolution_functions = deconvolution_functions))
saveRDS(bivariate_scenario, "/home/bncl_cb/rstudio/working/DeCovarT/simulations/results/bivariate_scenario.rds")

# easy_scenario <- DeCovarT::benchmark_deconvolution_algorithms(proportions=list("balanced"=c(0.5, 0.5)), n = 200,
#                                                     signature_matrices=matrix(c(20, 40, 40, 20), nrow = 2),
#                                                     deconvolution_functions = list("LBFGS"=list(FUN=DeCovarT::deconvolute_ratios_LBFGS),
#                                                                                    "lsei"=list(FUN=DeCovarT::deconvolute_ratios_deconRNASeq),
#                                                                                    "optim"=list(FUN=DeCovarT::deconvolute_ratios_basic_optim),
#                                                                                    "barrier"=list(FUN=DeCovarT::deconvolute_ratios_constrOptim),
#                                                                                    "DeCoVarT"=list(FUN=DeCovarT::deconvolute_ratios_DeCoVarT),
#                                                                                    "NLM"=list(FUN=DeCovarT::deconvolute_ratios_nlm),
#                                                                                    "SA"=list(FUN=DeCovarT::deconvolute_ratios_simulated_annealing)))
# saveRDS(easy_scenario, "/home/bncl_cb/rstudio/working/DeCovarT/simulations/results/easy_scenario.rds")
# 
# complex_scenario <- DeCovarT::benchmark_deconvolution_algorithms(proportions=list("highly_unbalanced"=c(0.95, 0.05)), n = 200,
#                                                        signature_matrices=matrix(c(20, 40, 40, 20), nrow = 2),
#                                                        deconvolution_functions = list("lsei"=list(FUN=DeCovarT::deconvolute_ratios_deconRNASeq),
#                                                                                       "optim"=list(FUN=DeCovarT::deconvolute_ratios_basic_optim),
#                                                                                       "barrier"=list(FUN=DeCovarT::deconvolute_ratios_constrOptim),
#                                                                                       "DeCoVarT"=list(FUN=DeCovarT::deconvolute_ratios_DeCoVarT),
#                                                                                       "NLM"=list(FUN=DeCovarT::deconvolute_ratios_nlm),
#                                                                                       "SA"=list(FUN=DeCovarT::deconvolute_ratios_simulated_annealing)))
# saveRDS(complex_scenario, "/home/bncl_cb/rstudio/working/DeCovarT/simulations/results/complex_scenario.rds")
# 


# RNGkind("L'Ecuyer-CMRG"); set.seed(3) # set an unique core
# small_overlap <- DeCovarT::benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.5, 0.5), "highly_unbalanced"=c(0.95, 0.05)), n = 500,
#                                                     signature_matrices=matrix(c(20, 40, 40, 20), nrow = 2), 
#                                                     diagonal_terms = list("homoscedasctic"= diag(c(1, 1))),
#                                                     deconvolution_functions = list("abbas"=list(FUN=DeCovarT::deconvolute_ratios_abbas),
#                                                                                    "lsei"=list(FUN=DeCovarT::deconvolute_ratios_monaco),
#                                                                                    "nlm"=list(FUN=DeCovarT::deconvolute_ratios_nlm),
#                                                                                    "DeCoVarT"=list(FUN=DeCovarT::deconvolute_ratios_DeCoVarT),
#                                                                                    "tricky optim"=list(FUN=DeCovarT::deconvolute_ratios_constrained_optim),
#                                                                                    "SA"=list(FUN=DeCovarT::deconvolute_ratios_simulated_annealing)))
# saveRDS(small_overlap, "/home/bncl_cb/rstudio/working/DeCovarT/simulations/results/small_overlap_version2.rds")
# 
# 
# 
# RNGkind("L'Ecuyer-CMRG"); set.seed(3) # set an unique core
# high_overlap <- DeCovarT::benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.5, 0.5), "highly_unbalanced"=c(0.95, 0.05)), n = 500,
#                                                    signature_matrices=matrix(c(20, 22, 22, 20), nrow = 2),
#                                                    diagonal_terms = list("homoscedasctic"= diag(c(1, 1))),
#                                                    deconvolution_functions = list("abbas"=list(FUN=DeCovarT::deconvolute_ratios_abbas),
#                                                                                   "lsei"=list(FUN=DeCovarT::deconvolute_ratios_monaco),
#                                                                                   "nlm"=list(FUN=DeCovarT::deconvolute_ratios_nlm),
#                                                                                   "DeCoVarT"=list(FUN=DeCovarT::deconvolute_ratios_DeCoVarT),
#                                                                                   "tricky optim"=list(FUN=DeCovarT::deconvolute_ratios_constrained_optim),
#                                                                                   "SA"=list(FUN=DeCovarT::deconvolute_ratios_simulated_annealing)))
# 
# saveRDS(high_overlap, "/home/bncl_cb/rstudio/working/DeCovarT/simulations/results/high_overlap_version2.rds")









