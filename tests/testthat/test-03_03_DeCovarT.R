test_that("Use numDeriv to check numerically derviative values", {
  set.seed(3)
  library(tensorA)
  mean_signature_matrix <- matrix(c(20, 40, 40, 20), nrow = 2); p <- c(0.5, 0.5)
  num_genes <- nrow(mean_signature_matrix); num_celltypes <- ncol(mean_signature_matrix)
  y <- mean_signature_matrix %*% p + rnorm(nrow(mean_signature_matrix)) # global gene expression, as linear combination
  Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
                 dim = c(num_genes,num_genes,num_celltypes))
  ##----------------------------------------------------------------
  ##                  test gradient log-likelihood                 -
  ##----------------------------------------------------------------
  jacobian_mapping_numerical <- numDeriv::grad(loglik_multivariate, p,
                                               method="Richardson", method.args=list(eps=1e-4, r=6),
                                               y=y, mean_signature_matrix=mean_signature_matrix, Sigma=Sigma) # additional arguments

  jacobian_mapping_theoretical <- gradient_loglik_unconstrained (p, y, mean_signature_matrix, Sigma)
  testthat::expect_equal(jacobian_mapping_numerical, jacobian_mapping_theoretical)
})


test_that("Check Maximum a posterior computation", {
  set.seed(3)
  # Define simulation parameters --------------------------------------------
  mean_signature_matrix <- matrix(c(20, 40, 40, 20), nrow = 2, 
                                  dimnames = list(paste0("gene_", 1:2), paste0("celltype_", 1:2)))
  p <- c(0.5, 0.5)
  num_genes <- nrow(mean_signature_matrix); num_celltypes <- ncol(mean_signature_matrix)
  Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
                 dim = c(num_genes,num_genes,num_celltypes), 
                 dimnames = list(paste0("gene_", 1:2), paste0("gene_", 1:2), paste0("celltype_", 1:2)))

  # Simulate accordingly ----------------------------------------------------
  simulated_data <- simulate_bulk_mixture ( signature_matrix = mean_signature_matrix,
                                            Sigma,
                                            p = c(0.5, 0.5),
                                            n = 1)
  y_simu <- simulated_data$Y; X_simu <- simulated_data$mean_signature_matrix
  rm(simulated_data)

# Compute MAP gaussian ----------------------------------------------------
MAP_two_ratios <- .map_gaussian_convolution (y = y_simu,
                                             mean_signature_matrix, 
                                             Sigma)
  
  
  testthat::expect_equal(jacobian_mapping_numerical, jacobian_mapping_theoretical)
})

test_that("compare the performance of several algorithms", {
  set.seed(3)
  library(dplyr)
  

# Define simulation parameters --------------------------------------------
  mean_signature_matrix <- matrix(c(20, 40, 40, 20), nrow = 2, 
                                  dimnames = list(paste0("gene_", 1:2), paste0("celltype_", 1:2)))
  p <- c(0.5, 0.5)
  num_genes <- nrow(mean_signature_matrix); num_celltypes <- ncol(mean_signature_matrix)
  Sigma <- array(c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
                 dim = c(num_genes,num_genes,num_celltypes), 
                 dimnames = list(paste0("gene_", 1:2), paste0("gene_", 1:2), paste0("celltype_", 1:2)))

# Simulate accordingly ----------------------------------------------------
  simulated_data <- simulate_bulk_mixture ( signature_matrix = mean_signature_matrix,
                                            Sigma,
                                            p = c(0.5, 0.5),
                                            n = 2000)
  y_simu <- simulated_data$Y; X_simu <- simulated_data$mean_signature_matrix

# Deconvolution of synthetic data -------------------
  inferred_ratios <- deconvolute_ratios_Newton_Raphson(
    y = y_simu[, 2, drop=FALSE],
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    true_ratios = c(0.5, 0.5),
    epsilon = 10^-4,
    itmax = 200
  )

  })
