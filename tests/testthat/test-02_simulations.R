##################################################################
##                   two-component simulation                   ##
##################################################################


test_that("small simulation testing", {
  simulation_two_genes <- simulate_bulk_mixture (signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2,
                                                                        dimnames = list(paste0("genes_", 1:2), paste0("cell_type_", 1:2))),
                                                           cov_tensor=array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2),
                                                                            dimnames = list(paste0("genes_", 1:2), paste0("genes_", 1:2), paste0("cell_type_", 1:2))), n=10)


})

test_that("test performance of several deconvolution algorithms", {
  two_genes_small_OVL_benchmark <- benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.5, 0.5)), n = 4,
                                                            signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2),
                                                            deconvolution_functions = list("lm"=list(FUN=deconvolute_ratios_abbas)))

  simulation_two_genes <- simulate_bulk_mixture (signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2,
                                                                         dimnames = list(paste0("genes_", 1:2), paste0("cell_type_", 1:2))),
                                                 cov_tensor=array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2),
                                                                  dimnames = list(paste0("genes_", 1:2), paste0("genes_", 1:2), paste0("cell_type_", 1:2))), n=1)
  
  

 
  estimated_p_DeCoVarT <- deconvolute_ratios_DeCoVarT(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                             Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_optim <- deconvolute_ratios_basic_optim(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                      Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_tricky_optim <- deconvolute_ratios_constrained_optim(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                      Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))

  estimated_p_constr<- deconvolute_ratios_constrOptim(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_LBFGS <- deconvolute_ratios_LBFGS(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                      Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))

  estimated_p_nlm <- deconvolute_ratios_nlm(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_SA <- deconvolute_ratios_simulated_annealing(y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                            Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))



  estimated_p_abbas <- deconvolute_ratios_abbas  (y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                 Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_deconRNASeq <- deconvolute_ratios_deconRNASeq  (y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                  Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_monaco <- deconvolute_ratios_monaco  (y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                  Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))
  estimated_p_nnls <- deconvolute_ratios_nnls  (y = simulation_two_genes$Y, X=simulation_two_genes$X[,,1] %>% as.matrix(),
                                                  Sigma = array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2), dim = c(2,2,2)), true_ratios = c(0.5, 0.5))

})


test_that("understand mistakes", {
  erreur_output <- readRDS("./simulations/erreurs/erreur_1_function_DeCoVarT.rds")
  wrong_estimate <- do.call(deconvolute_ratios_DeCoVarT, erreur_output)
})

test_that("benchmark deconvolution", {
  RNGkind("L'Ecuyer-CMRG"); set.seed(3) # set an unique core
  complex_scenario <- benchmark_deconvolution_algorithms_two_genes(proportions=list("unbalanced"=c(0.95, 0.05)), n = 4, 
                                                                corr_sequence = c(-0.8,  0.8),
                                                                signature_matrix=matrix(c(20, 22, 22, 20), nrow = 2),
                                                                deconvolution_functions = list("DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT)))
  
  
  saveRDS(easy_scenario, "./simulations/results/toy_example_easy_scenario.rds")

  complex_scenario <- benchmark_deconvolution_algorithms_two_genes(proportions=list("highly_unbalanced"=c(0.95, 0.05)), n = 500,
                                                      signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2),
                                                      deconvolution_functions = list("lsei"=list(FUN=deconvolute_ratios_deconRNASeq),
                                                                                     "optim"=list(FUN=deconvolute_ratios_basic_optim),
                                                                                     "barrier"=list(FUN=deconvolute_ratios_constrOptim),
                                                                                     "LBFGS"=list(FUN=deconvolute_ratios_LBFGS),
                                                                                     "DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT),
                                                                                     "NLM"=list(FUN=deconvolute_ratios_nlm),
                                                                                     "SA"=list(FUN=deconvolute_ratios_simulated_annealing)))
  saveRDS(complex_scenario, "./simulations/results/complex_scenario.rds")

  small_overlap <- benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.5, 0.5), "highly_unbalanced"=c(0.95, 0.05)), n = 500,
                                                      signature_matrix=matrix(c(20, 40, 40, 20), nrow = 2),
                                                      deconvolution_functions = list("abbas"=list(FUN=deconvolute_ratios_abbas),
                                                                                     "nnls"=list(FUN=deconvolute_ratios_nnls),
                                                                                     "abbas"=list(FUN=deconvolute_ratios_abbas),
                                                                                     "barrier"=list(FUN=deconvolute_ratios_constrOptim),
                                                                                     "DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT)))
  saveRDS(small_overlap, "./simulations/results/small_overlap.rds")


  high_overlap <- benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.5, 0.5), "highly_unbalanced"=c(0.95, 0.05)), n = 500,
                                                     signature_matrix=matrix(c(20, 24, 24, 20), nrow = 2),
                                                     deconvolution_functions = list("abbas"=list(FUN=deconvolute_ratios_abbas),
                                                                                    "nnls"=list(FUN=deconvolute_ratios_nnls),
                                                                                    "abbas"=list(FUN=deconvolute_ratios_abbas),
                                                                                    "barrier"=list(FUN=deconvolute_ratios_constrOptim),
                                                                                    "DeCoVarT"=list(FUN=deconvolute_ratios_DeCoVarT)))
  saveRDS(high_overlap, "./simulations/results/high_overlap.rds")

})














# saveRDS(two_genes_small_OVL, file = "../simulations/two_genes_small_OVL.rds")

metric_plots <- purrr::map(two_genes_small_OVL, "metric_plots") %>%
  gridExtra::marrangeGrob(ncol=1, nrow=1, top = "")

ggsave("../results/ComplexHeatmap/two_genes_small_OVL_multivariate.pdf",
       metric_plots, width = 12, height = 14,dpi = 600)

dp_plots <- two_genes_small_OVL$balanced$simulations_dp

dp_plots_small_OVL <- cowplot::plot_grid(
  plotlist = list(dp_plots %>%
                    filter(correlation_celltype1==-0.8 & correlation_celltype2 == -0.8 & variance=="heteroscedastic") %>%
                    pull(density_plot) %>% purrr::pluck(1),
                  dp_plots %>%
                    filter(correlation_celltype1==-0.8 & correlation_celltype2 == 0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0.8 & correlation_celltype2 == 0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0 & correlation_celltype2 == 0 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1)) %>% purrr::simplify(), align = "hv", axis = "tblr", nrow = 2)
ggsave("../results/density_plots/dp_small_OVL_balanced_two_genes.pdf", dp_plots_small_OVL,
       dpi = 600, width = 14, height = 14)



# simulation with high OVL
two_genes_heteroscedastic_high_OVL <- benchmark_deconvolution_algorithms_two_genes(proportions=list("balanced"=c(0.5, 0.5), "small unbalanced"=c(0.75, 0.25),"highly unbalanced"=c(0.95, 0.05)), n = 2000,
                                                                         signature_matrix=matrix(c(20, 22, 22, 20), nrow = 2),
                                                                         deconvolution_functions = list("lm" = list(FUN=deconvolute_ratios_abbas,additionnal_parameters=NULL)))

saveRDS(two_genes_heteroscedastic_high_OVL, file = "../simulations/two_genes_heteroscedastic_high_OVL.rds")
metric_plots <- purrr::map(two_genes_heteroscedastic_high_OVL, "metric_plots") %>%
  gridExtra::marrangeGrob(ncol=1, nrow=1, top = "")
ggsave("../results/ComplexHeatmap/two_genes_high_OVL.pdf",
       metric_plots, width = 14, height = 18,dpi = 600)

dp_plots <- two_genes_heteroscedastic_high_OVL$balanced$simulations_dp
dp_plots_high_OVL <- cowplot::plot_grid(
  plotlist = list(dp_plots %>%
                    filter(correlation_celltype1==-0.8 & correlation_celltype2 == -0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% purrr::pluck(1),
                  dp_plots %>%
                    filter(correlation_celltype1==-0.8 & correlation_celltype2 == 0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0.8 & correlation_celltype2 == 0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0.8 & correlation_celltype2 == -0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==-0.4 & correlation_celltype2 == -0.4 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0 & correlation_celltype2 == 0 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1)) %>% purrr::simplify(), align = "hv", axis = "tblr", ncol = 2)
ggsave("../results/density_plots/dp_high_OVL_balanced_two_genes.pdf", dp_plots_high_OVL,
       dpi = 600, width = 14, height = 18)

dp_plots <- two_genes_heteroscedastic_high_OVL$`highly unbalanced`$simulations_dp
dp_plots_high_OVL_highly_unbalanced <- cowplot::plot_grid(
  plotlist = list(dp_plots %>%
                    filter(correlation_celltype1==-0.8 & correlation_celltype2 == -0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% purrr::pluck(1),
                  dp_plots %>%
                    filter(correlation_celltype1==-0.8 & correlation_celltype2 == 0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0.8 & correlation_celltype2 == 0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0.8 & correlation_celltype2 == -0.8 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==-0.4 & correlation_celltype2 == -0.4 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1),
                  dp_plots %>%
                    filter(correlation_celltype1==0 & correlation_celltype2 == 0 & variance=="homoscedastic") %>%
                    pull(density_plot) %>% magrittr::extract2(1)) %>% purrr::simplify(), align = "hv", axis = "tblr", ncol = 2)
ggsave("../results/density_plots/dp_plots_high_OVL_highly_unbalanced.pdf", dp_plots_high_OVL_highly_unbalanced,
       dpi = 600, width = 14, height = 18)

mean_sd <- list(mean = ~mean(.x, na.rm = TRUE),
                sd = ~sd(.x, na.rm = TRUE))

simulation_metrics_summary <- two_genes_heteroscedastic_high_OVL %>% purrr::map_dfr("simulation_metrics")  %>%
  group_by(proportions, correlation_celltype1, correlation_celltype2, variance, deconvolution_name) %>%
  summarise(across(c(model_mse, model_rmse, model_coef_determination, model_mae, celltype_1, celltype_2), mean_sd)) %>%
  arrange(correlation_celltype1, correlation_celltype2, desc(proportions), desc(variance))



# to save the distributions in pdf objects
# upper_part <- cowplot::plot_grid(density_plots$individual_plots$UR_0_skewness_0_OVL_0.16 + labs(tag="A") +
#                                    theme(plot.tag = element_text(size=18, face = "bold")),
#                                  two_components_balanced_overlapping_microbenchmark_computations$init_focus_time_plots$kmeans +
#                                    theme(strip.text = element_blank(), plot.title = element_blank(), plot.tag = element_text(size=18, face = "bold"),
#                                          legend.text = element_text(size=16), axis.title=element_text(size=16)) +
#                                    labs(tag = "B", x="Number of observations (log10)", y="Time in seconds (log10)"),
#                                  ncol = 2, align = 'h', axis="tblr")
#
# two_components_balanced_overlapping <- gridExtra::arrangeGrob(upper_part,
#                                                               plot_boxplots_parameters(benchmark_scores_2000_observations_UR_0_skewness_0_OVL_0.16_prop_outliers_0$plots$data, p=c(0.5, 0.5), mu=c(0, 4), sigma=c(2,2)) +
#                                                                 theme(plot.tag = element_text(size=18, face = "bold")) + labs(tag = "C"),
#                                                               top="", nrow = 2, heights = c(1.2, 2), padding = unit(1, "line"))
# ggsave("paper_figures/two_components_balanced_overlapping.pdf", two_components_balanced_overlapping, width = 15, height = 18,dpi = 600)





##################################################################
##                  multi-component simulation                  ##
##################################################################


# library(bnlearn)
# data(gaussian.test)
# dag = model2network("[A][B][E][G][C|A:B][D|B][F|A:D:E:G]")
# bn = bn.fit(dag, gaussian.test)
# mvn = gbn2mvnorm(bn)
# bn2 = mvnorm2gbn(dag, mu = mvn$mu, sigma = mvn$sigma)
# all.equal(bn, bn2)



############################################################################
############################################################################
###                                                                      ###
###                     TEST NUMERICALLY DERIVATIVES                     ###
###                                                                      ###
############################################################################
############################################################################



#################################################################
##                  first part of the Hessian                  ##
#################################################################
library(tensorA) # load package
J_log <- to.tensor(1:3,c("i"=3)) # Jacobian of the original log-likelihood
H_psi <- to.tensor(1:12,c(j=2,k=2,"i"=3)) # Hessian of the mapping function

# Einstein summation
hess_einstein <- J_log %e% H_psi # %e% refers to the Einstein operator, repeated indexes (here i)
# in both tensors are coupled and summed over, of note, any number of tensors can be passed
# %e% operator, provided you respect the shape configuration

# multivariate tensor product
hess_tensor_product <- mul.tensor(J_log,"i",H_psi,"i") # the standard tensor product,
# in which you have to indicate explicitly the indexes to aggregate

#  with the tensor package
library(tensor)
# in opposition to the mult product, you instead indicate the distinct indexes
hess_tensor <- tensor(A = J_log, B = H_psi, alongA = 1, alongB = 3)

# with base R operators
hess_base <- matrix(0, nrow = 2, ncol = 2, dimnames = )
for (i in 1:3) {
  hess_base <- hess_base +
    (J_log[i] %>% as.numeric()) * (H_psi[,,i] %>% to.matrix.tensor(j="I2"))
}

# check that all these elements are equal
testthat::expect_equal(hess_rieman, hess_tensor %>% to.tensor()) # pass
testthat::expect_equal(hess_rieman, hess_tensor_product) # pass
# pass but requires cumbersome operation
testthat::expect_equal(hess_rieman, hess_base %>%
                         to.tensor() %>% extract2(I1 =~"j", I2 =~"k"))

##################################################################
##                  second part of the Hessian                  ##
##################################################################

J_psi <- to.tensor(1:6,c("i"=3, "j"=2))
H_log <- to.tensor(1:9, c("i"=3, "l"=3))

# with base R
hess_base_R <- t(J_psi) %*% H_log %*% J_psi
# %*% designs here the dot product

# with the Einstein summation (needs an index change)
hess_einstein_second <- J_psi %e% H_log %e% J_psi[[i =~"l", j =~"k"]]
testthat::expect_equal(hess_base_R,
                       hess_einstein_second %>% matrix(nrow=2, ncol=2)) # pass



# testthat::expect_equal(J_psi %e% H_log %e% J_psi,
# mul.tensor(J_psi, i="i", H_log, j="i") %>% mul.tensor(i="k", J_psi, j="k"))
#
#
# J_psi <- to.tensor(1:6,c("i"=3, "k"=2))
# H_log <- to.tensor(1:9, c("i"=3, "j"=3))
# riemann.tensor(J_psi, H_log[[i =~"^i"]], J_psi[[k =~"^k"]])


##################################################################
##                  test numerical derivatives                  ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")
equibalanced_ratios <- c(1/3, 1/3, 1/3)
testthat::expect_equal(equibalanced_ratios,
                       equibalanced_ratios %>% inverse_mapping_function() %>% mapping_function())

# check Jacobian of the mapping function
jacobian_mapping_theoretical <- jacobian_mapping_function(c(0, 0, 0))
jacobian_mapping_numerical <- numDeriv::jacobian(mapping_function, c(0, 0, 0), method="simple", method.args=list(eps=1e-6))
testthat::expect_equal(jacobian_mapping_theoretical, jacobian_mapping_numerical, tolerance = 1e-6) # pass



jacobian_mapping_theoretical <- jacobian_mapping_function(inverse_mapping_function(c(0.90, 0.05, 0.05)))
# privilege the Ridchardon and its two-way estimation of derivative #
jacobian_mapping_numerical <- numDeriv::jacobian(mapping_function, inverse_mapping_function(c(0.90, 0.05, 0.05)),
                                                 method="Richardson", method.args=list(eps=1e-6))
testthat::expect_equal(jacobian_mapping_theoretical, jacobian_mapping_numerical, tolerance = 1e-6) # pass


##----------------------------------------------------------------
##                  test gradient log-likelihood                 -
##----------------------------------------------------------------

set.seed(3)
X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
               dim = c(num_genes,num_genes,num_celltypes))


jacobian_mapping_numerical <- numDeriv::grad(loglik_multivariate, p,
                                             method="Richardson", method.args=list(eps=1e-4, r=6),
                                             y=y, X=X, Sigma=Sigma) # additional arguments

jacobian_mapping_theoretical <- gradient_loglik_unconstrained (p, y, X, Sigma)
testthat::expect_equal(jacobian_mapping_numerical, jacobian_mapping_theoretical)

#################################################################
##               test constrained log-likelihood               ##
#################################################################
set.seed(3); library(tensorA)
X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
               dim = c(num_genes,num_genes,num_celltypes))
theta <- inverse_mapping_function(p)

Sigma_global <- .compute_global_variance(p, Sigma)

theta_vect <- seq(-30, 30, 0.5)
log_vect <- purrr::map_dbl(theta_vect, loglik_multivariate_constrained, y=y, X=X, Sigma=Sigma)
plot(log_vect ~theta_vect)

loglik_multivariate_constrained(theta_vect[1], )

# matrix product of the gradient unconstrained evaluated in p times Jacobian of p
grad_constrained_mapping_numerical <- numDeriv::grad(loglik_multivariate_constrained, theta,
                                                     method="Richardson", method.args=list(eps=1e-4, r=6),
                                                     y=y, X=X, Sigma=Sigma) # additional arguments
grad_constrained_mapping_theoretical <- gradient_loglik_constrained (theta, y, X, Sigma)
testthat::expect_equal(grad_constrained_mapping_numerical, grad_constrained_mapping_theoretical %>% as.numeric())

###  more complex derivative
set.seed(3)
X <- matrix(c(20, 40, 50, 10, 40, 20, 60, 40, 25), nrow = 3); p <- c(0.2, 0.3, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
# y <- X %*% p # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.4, 0.8, 1, 0.2, 0.4, 0.2, 1,
                 2, -0.2, -0.2, -0.2, 2, -0.2, -0.2, -0.2, 2,
                 1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1),
               dim = c(num_genes,num_genes,num_celltypes))
theta <- inverse_mapping_function(p)

# matrix product of the gradient unconstrained evaluated in p times Jacobian of p
grad_constrained_mapping_numerical <- numDeriv::grad(loglik_multivariate_constrained, theta,
                                                     method="Richardson", method.args=list(eps=1e-12, r=4),
                                                     y=y, X=X, Sigma=Sigma) # additional arguments
grad_constrained_mapping_theoretical <- gradient_loglik_constrained (theta, y, X, Sigma)
testthat::expect_equal(grad_constrained_mapping_numerical, c(grad_constrained_mapping_theoretical))

##################################################################
##                 loglik_hessian_unconstrained                 ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")

set.seed(3)
X <- matrix(c(20, 40, 50, 10, 40, 20, 60, 40, 25), nrow = 3); p <- c(0.2, 0.3, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X))# global gene expression, as linear combination
# y <- X %*% p # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.4, 0.8, 1, 0.2, 0.4, 0.2, 1,
                 2, -0.2, -0.2, -0.2, 2, -0.2, -0.2, -0.2, 2,
                 1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1),
               dim = c(num_genes,num_genes,num_celltypes))

hessian_loglik_numerical <- numDeriv::hessian(loglik_multivariate, p,
                                              method="Richardson", method.args=list(eps=1e-12, r=4),
                                              y=y, X=X, Sigma=Sigma) # additional arguments
# jacobian_grad_numerical <- numDeriv::jacobian(gradient_loglik_unconstrained, p,
#                                               method="Richardson", method.args=list(eps=1e-12, r=4),
#                                               y=y, X=X, Sigma=Sigma)
# testthat::expect_equal(hessian_loglik_numerical, jacobian_grad_numerical)


hessian_loglik_theoretical <- hessian_loglik_unconstrained (p, y, X, Sigma)


set.seed(3)
X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
num_genes <- nrow(X); num_celltypes <- ncol(X)
y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
# y <- X %*% p
Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
               dim = c(num_genes,num_genes,num_celltypes))

hessian_loglik_theoretical <- hessian_loglik_unconstrained (p, y, X, Sigma)
numDeriv::hessian(loglik_multivariate, p,
                  method="Richardson", method.args=list(eps=1e-12, r=4),
                  y=y, X=X, Sigma=Sigma)

##################################################################
##                 hessian mapping                 ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")
p <- c(0.2, 0.3, 0.4, 0.1); theta <- inverse_mapping_function(p)

# check Hessian of the mapping function
hessian_mapping_theoretical <- hessian_mapping_function(theta)
hessian_mapping_numerical <- array(0, dim = c(length(theta), length(theta), length(p)))
for (i in 1:length(p)) {
  hessian_mapping_numerical[,,i] <- numDeriv::hessian(mapping_function_univariate, theta,
                                                      method="Richardson", method.args=list(eps=1e-6), index=i)
}
testthat::expect_equal(hessian_mapping_theoretical, hessian_mapping_numerical, tolerance = 1e-6) # pass



##################################################################
##                 loglik_hessian_constrained                 ##
##################################################################
library(dplyr); source("R/utils.R")
source("R/03_01_deconvolution_algorithms.R"); source("R/03_01_bis_deconvolution_algorithms.R")

set.seed(3)
X <- matrix(c(20, 40, 50, 10, 40, 20, 60, 40, 25), nrow = 3); p <- c(0.2, 0.3, 0.5)
theta <- inverse_mapping_function(p)
num_genes <- nrow(X); num_celltypes <- ncol(X)
# y <- X %*% p + rnorm(nrow(X))# global gene expression, as linear combination
y <- X %*% p # global gene expression, as linear combination
Sigma <- array(c(1, 0.8, 0.4, 0.8, 1, 0.2, 0.4, 0.2, 1,
                 2, -0.2, -0.2, -0.2, 2, -0.2, -0.2, -0.2, 2,
                 1, 0.4, 0.4, 0.4, 1, 0.4, 0.4, 0.4, 1),
               dim = c(num_genes,num_genes,num_celltypes))



hessian_constrained_mapping_numerical <- numDeriv::hessian(loglik_multivariate_constrained, theta,
                                                           method="Richardson", method.args=list(eps=1e-12, r=4),
                                                           y=y, X=X, Sigma=Sigma) # additional arguments
hessian_constrained_mapping_theoretical <- hessian_loglik_constrained (theta, y, X, Sigma)
testthat::expect_equal(hessian_constrained_mapping_numerical, hessian_constrained_mapping_theoretical)






