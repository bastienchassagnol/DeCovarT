test_that("Check Maximum a Posteriori computation", {
  # Define simulation parameters --------------------------------------------
  mean_signature_matrix <- matrix(
    c(20, 40, 40, 20),
    nrow = 2,
    dimnames = list(paste0("gene_", 1:2), paste0("celltype_", 1:2))
  )
  p <- c(0.5, 0.5)
  num_genes <- nrow(mean_signature_matrix)
  num_celltypes <- ncol(mean_signature_matrix)
  Sigma <- array(
    c(1, 0.8, 0.8, 1, 2, -0.2, -0.2, 2),
    dim = c(num_genes, num_genes, num_celltypes),
    dimnames = list(
      paste0("gene_", 1:2),
      paste0("gene_", 1:2),
      paste0("celltype_", 1:2)
    )
  )

  # Simulate accordingly ----------------------------------------------------
  simulated_data <- withr::with_seed(
    3L,
    simulate_bulk_mixture(
      signature_matrix = mean_signature_matrix,
      Sigma,
      p = c(0.5, 0.5),
      n = 1
    )
  )
  y_simu <- simulated_data$Y
  X_simu <- simulated_data$mean_signature_matrix
  rm(simulated_data)

  # Compute MAP gaussian ----------------------------------------------------
  MAP_two_ratios <- .map_gaussian_convolution(
    y = y_simu,
    mean_signature_matrix,
    Sigma
  )

  testthat::expect_equal(
    jacobian_mapping_numerical,
    jacobian_mapping_theoretical
  )
})
