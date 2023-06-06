# test_that("assert numerically tensor product", {
#   #################################################################
#   ##                  first part of the Hessian                  ##
#   #################################################################
#   library(tensorA) # load package
#   J_log <- to.tensor(1:3,c("i"=3)) # Jacobian of the original log-likelihood
#   H_psi <- to.tensor(1:12,c(j=2,k=2,"i"=3)) # Hessian of the mapping function
#
#   # Einstein summation
#   hess_einstein <- J_log %e% H_psi # %e% refers to the Einstein operator, repeated indexes (here i)
#   # in both tensors are coupled and summed over, of note, any number of tensors can be passed
#   # %e% operator, provided you respect the shape configuration
#
#   # multivariate tensor product
#   hess_tensor_product <- mul.tensor(J_log,"i",H_psi,"i") # the standard tensor product,
#   # in which you have to indicate explicitly the indexes to aggregate
#
#   #  with the tensor package
#   library(tensor)
#   # in opposition to the mult product, you instead indicate the distinct indexes
#   hess_tensor <- tensor(A = J_log, B = H_psi, alongA = 1, alongB = 3)
#
#   # with base R operators
#   hess_base <- matrix(0, nrow = 2, ncol = 2, dimnames = )
#   for (i in 1:3) {
#     hess_base <- hess_base +
#       (J_log[i] %>% as.numeric()) * (H_psi[,,i] %>% to.matrix.tensor(j="I2"))
#   }
#
#
#
#   # check that all these elements are equal
#   testthat::expect_equal(hess_rieman, hess_tensor %>% to.tensor()) # pass
#   testthat::expect_equal(hess_rieman, hess_tensor_product) # pass
#   # pass but requires cumbersome operation
#   testthat::expect_equal(hess_rieman, hess_base %>%
#                            to.tensor() %>% extract2(I1 =~"j", I2 =~"k"))
#
#   ##################################################################
#   ##                  second part of the Hessian                  ##
#   ##################################################################
#
#   J_psi <- to.tensor(1:6,c("i"=3, "j"=2))
#   H_log <- to.tensor(1:9, c("i"=3, "l"=3))
#
#   # with base R
#   hess_base_R <- t(J_psi) %*% H_log %*% J_psi
#   # %*% designs here the dot product
#
#   # with the Einstein summation (needs an index change)
#   hess_einstein_second <- J_psi %e% H_log %e% J_psi[[i =~"l", j =~"k"]]
#   testthat::expect_equal(hess_base_R,
#                          hess_einstein_second %>% matrix(nrow=2, ncol=2)) # pass
# })
#
#
# test_that("assert numerically gradient and Hessian of the reparametrised function", {
#   set.seed(3); library(tensorA)
#   X <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
#   num_genes <- nrow(X); num_celltypes <- ncol(X)
#   y <- X %*% p + rnorm(nrow(X)) # global gene expression, as linear combination
#   Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
#                  dim = c(num_genes,num_genes,num_celltypes))
#   theta <- inverse_mapping_function(p)
#   ##----------------------------------------------------------------
#   ##                  test gradient log-likelihood                 -
#   ##----------------------------------------------------------------
#   jacobian_mapping_numerical <- numDeriv::grad(loglik_multivariate, p,
#                                                method="Richardson", method.args=list(eps=1e-4, r=6),
#                                                y=y, X=X, Sigma=Sigma) # additional arguments
#
#   jacobian_mapping_theoretical <- gradient_loglik_unconstrained (p, y, X, Sigma)
#   testthat::expect_equal(jacobian_mapping_numerical, jacobian_mapping_theoretical)
#
#   #################################################################
#   ##               test constrained log-likelihood               ##
#   #################################################################
#   # matrix product of the gradient unconstrained evaluated in p times Jacobian of p
#   grad_constrained_mapping_numerical <- numDeriv::grad(loglik_multivariate_constrained, theta,
#                                                        method="Richardson", method.args=list(eps=1e-4, r=6),
#                                                        y=y, X=X, Sigma=Sigma) # additional arguments
#   grad_constrained_mapping_theoretical <- gradient_loglik_constrained (theta, y, X, Sigma)
#   testthat::expect_equal(grad_constrained_mapping_numerical, grad_constrained_mapping_theoretical %>% as.numeric())
#
# ##################################################################
# ##                 loglik_hessian_unconstrained                 ##
# ##################################################################
# hessian_loglik_numerical <- numDeriv::hessian(loglik_multivariate, p,
#                                               method="Richardson", method.args=list(eps=1e-12, r=4),
#                                               y=y, X=X, Sigma=Sigma) # additional arguments
# jacobian_grad_numerical <- numDeriv::jacobian(gradient_loglik_unconstrained, p,
#                                               method="Richardson", method.args=list(eps=1e-12, r=4),
#                                               y=y, X=X, Sigma=Sigma)
# testthat::expect_equal(hessian_loglik_numerical, jacobian_grad_numerical)
# hessian_loglik_theoretical <- hessian_loglik_unconstrained (p, y, X, Sigma)
#
# ##################################################################
# ##                 loglik_hessian_constrained                 ##
# ##################################################################
# hessian_constrained_mapping_numerical <- numDeriv::hessian(loglik_multivariate_constrained, theta,
#                                                            method="Richardson", method.args=list(eps=1e-12, r=4),
#                                                            y=y, X=X, Sigma=Sigma) # additional arguments
# hessian_constrained_mapping_theoretical <- hessian_loglik_constrained (theta, y, X, Sigma)
# testthat::expect_equal(hessian_constrained_mapping_numerical, hessian_constrained_mapping_theoretical)
# })
#
#
# test_that("compare the performance of several algorithms", {
#   set.seed(3); library(dplyr)
#   X <- matrix(c(20, 40, 40, 20), nrow = 2, dimnames = list(paste0("gene_", 1:2), paste0("celltype_", 1:2)))
#   p <- c(0.5, 0.5)
#   num_genes <- nrow(X); num_celltypes <- ncol(X)
#   y <- X %*% p + rnorm(nrow(X)); y <- as.vector(y) # global gene expression, as linear combination
#   Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
#                  dim = c(num_genes,num_genes,num_celltypes))
#   theta <- inverse_mapping_function(p)
#
#   ###  compare the estimation of several estimation methods
#   basic_estimates <- deconvolute_ratios_basic_optim (y, X, Sigma)
#   nnls_estimates <- deconvolute_ratios_nnls(y, X)
#
#
#
#   ###  understand recurrent mistakes
#   decovart_error_estimates <- do.call(deconvolute_ratios_DeCoVarT, erreur_1_function_DeCoVarT)
#
#
#
#   # high dimensional simulations
#   set.seed(3)
#   highly_overlapping_parameter <- MixSim::MixSim(BarOmega = 0.01,
#                                                  K=2, p=20, sph = FALSE, hom = FALSE,
#                                                  ecc = 0.9, PiLow = 0.05, int = c(3, 50))
#   highly_overlapping_parameter_formatted <- list(p=highly_overlapping_parameter$Pi,
#                                                  mu=t(highly_overlapping_parameter$Mu),
#                                                  sigma=highly_overlapping_parameter$S) %>%
#     enforce_parameter_identifiability()
#
#   X <- highly_overlapping_parameter_formatted$mu;
#   colnames(X) <- paste0("celltype_", 1:ncol(X)); row.names(X) <- paste0("gene_", 1:nrow(X))
#   # p <- c(0.05, 0.9, 0.05)
#   Sigma <- highly_overlapping_parameter_formatted$sigma; dimnames(Sigma) <- list(NULL, NULL, paste0("celltype_", 1:ncol(X)))
#   p <- highly_overlapping_parameter_formatted$p
#   y <- X %*% p + rnorm(nrow(X)); y <- as.vector(y) # global gene expression, as linear combination
#   theta <- inverse_mapping_function(p)
#
#   nnls_estimates <- deconvolute_ratios_nnls(y, X)
#   decovart_error_estimates <- deconvolute_ratios_DeCoVarT (y, X, Sigma)
#   deconvolute_ratios_constrOptim(y, X, Sigma)
#
#
#
#   simulated_data <- simulate_bulk_mixture (X, Sigma, proportions=p, n=10^4)
#   global_mu <- p[1]*X[,1] + p[2]*X[,2]
#   global_sigma <- p[1]^2*Sigma[,,1] + p[2]^2*Sigma[,,2]
#
#   y_simu <- simulated_data$Y; X_simu <- simulated_data$X
#
#   # test_metrics_Decovart <- purrr::map_dfr(y_simu[,1:10] %>% dplyr::as_tibble(),
#   #                                         deconvolute_ratios_DeCoVarT, X, Sigma)
#   test_metrics_nlm <- purrr::map_dfr(y_simu[,1:10] %>% dplyr::as_tibble(),
#                                           deconvolute_ratios_nlm, X, Sigma)
#   test_metrics_lsei <- purrr::map_dfr(y_simu[,1:10] %>% dplyr::as_tibble(),
#                                           deconvolute_ratios_deconRNASeq, X)
#
#
#   mapping_tested <- seq(-4, 4, by = 0.01)
#   log_lik_vect <- purrr::map_dbl(mapping_tested, loglik_multivariate_constrained, y_simu[,4], X, Sigma)
#   plot(log_lik_vect ~ mapping_tested, pch=1)
#
#   })
