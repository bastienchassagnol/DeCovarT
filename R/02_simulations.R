#' Simulate n bulk mixtures, knowing p, mu and cov individual variances
#'
#' @param signature_matrix
#' @param cov_tensor
#' @param proportions respective proportions of cell types
#' @param n number of samples
#'
#' @import ggplot2
#' @return list with X: simulated purified expression profiles (as an array) and Y: the reconstituted bulk mixtures
#' @export


simulate_bulk_mixture <- function(signature_matrix, cov_tensor,
                                  proportions=rep(1/ncol(signature_matrix), ncol(signature_matrix)), n=500) {

  ##################################################################
  ##                        check validity                        ##
  ##################################################################
  if (!all.equal(row.names(signature_matrix), union(dimnames(cov_tensor)[[1]], dimnames(cov_tensor)[[2]])))
    stop("Some of the genes are distinct between expected and covariance expression")
  else if (!all.equal(colnames(signature_matrix), dimnames(cov_tensor)[[3]]))
    stop("Cell types differ between expected and covariance expression")

  ##################################################################
  ##            remove undefined covariance structures            ##
  ##################################################################
  # first remove all cell types associated to a undefinite cov matrix
  valid_celltypes <- purrr::discard(colnames(signature_matrix), ~ cov_tensor %>% magrittr::extract(,,.x) %>% is.na() %>% any())
  # second, remove all cell types associated to a non-positive definite matrix
  valid_celltypes <- purrr::keep(valid_celltypes, ~ cov_tensor %>% magrittr::extract(,,.x) %>% is_positive_definite())
  # only keep cell_types associated to sound covariance structures
  cov_tensor <- cov_tensor[,,valid_celltypes]; signature_matrix <- signature_matrix[, valid_celltypes]

  ##### simulation
  valid_celltypes <- colnames(signature_matrix)
  names_genes <- row.names(signature_matrix)
  Y <- matrix(0, nrow = nrow(signature_matrix), ncol = n, dimnames = list(names_genes, paste0("sample_", 1:n)))
  X <- array(0, c(nrow(signature_matrix), ncol(signature_matrix),n),
             dimnames = list(names_genes, valid_celltypes, colnames(Y))) # store cell-type specific expression

  for (cell_name in valid_celltypes) { # generate individually each cell distribution
    mean_parameter <- signature_matrix[, cell_name]; covariance_parameter <- cov_tensor [,,cell_name]
    expression_per_celltype <- MASS::mvrnorm(n = n, mu=mean_parameter, Sigma=covariance_parameter,
                                             tol = 1e-12, empirical = FALSE) # simulate a draw from the covariance matrix
    X[,cell_name, ] <- t(expression_per_celltype)
  }
  
  # Y_test <- matrix(0, nrow = nrow(signature_matrix), ncol = n, dimnames = list(names_genes, paste0("sample_", 1:n)))
  # for (i in 1:n) {
  #   Y_test[,i] <- X[,,i] %*% proportions
  # }
  
  Y <- tensor::tensor(proportions, B = X, alongA = 1, alongB = 2) # tensor product does directly the computation X %*% p
  return(list(X=X, Y=Y))
}




#' Benchmark two-component cell populations
#'
#' @param proportions
#' @param signature_matrices
#' @param n
#' @param score_variable
#' @param scaled
#' @param deconvolution_functions
#'
#' @return
#' @export


benchmark_deconvolution_algorithms_two_genes <- function(proportions=list("balanced"=c(0.5, 0.5), "small unbalanced"=c(0.6, 0.4),"highly unbalanced"=c(0.05, 0.95)),
                                                         signature_matrices=list("small OVL" = matrix(c(20, 40, 40, 20), nrow = 2)), corr_sequence = seq(-0.8, 0.8, 0.2),
                                                         diagonal_terms = list("homoscedasctic"= diag(c(1, 1)), "heteroscedastic" = diag(c(1, 2) %>% sqrt())), 
                                                         deconvolution_functions=list("lm" = list(FUN=deconvolute_ratios_abbas,additionnal_parameters=NULL)),
                                                         n=200, scaled=FALSE) {
  ##################################################################
  ##            iterate over covariance structures            ##
  ##################################################################
  num_celltypes <- ncol(signature_matrices[[1]]); num_genes <- nrow(signature_matrices[[1]])
  signature_matrices <- purrr::map(signature_matrices, function (.x) {
    dimnames(.x) <- list(paste0("gene_", 1:num_genes), paste0("celltype_", 1:num_celltypes))
    return(.x)
  })
  id_scenario <- 1 # initiate ID scenario
  simulations <- purrr::imap_dfr(signature_matrices, function(mu, name_distance_centroids) {
    simulations_per_ratios <- purrr::imap_dfr(proportions, function(p, name_balance) {
      simulation_metrics <- tibble::tibble()
      for (corr_celltype1 in corr_sequence) {
        message(paste0("We are at equilibrium scenario: ", name_balance, " with correlation for celltype 1: ", corr_celltype1,
                    " and distance to centroids: ", name_distance_centroids, "."))
        for (corr_celltype2 in corr_sequence) {
          simulation_metrics <- simulation_metrics %>% dplyr::bind_rows(
            purrr::imap_dfr(diagonal_terms, function(.diag, .name) {
              ##------------------------------
              ##  generate covariance matrix
              ##------------------------------
              corr_matrix <- array(0, dim = c(num_genes,num_genes,num_celltypes),
                                   dimnames = list(paste0("gene_", 1:num_genes),
                                                   paste0("gene_", 1:num_genes),
                                                   paste0("celltype_", 1:num_celltypes)))
              cov_tensor <- corr_matrix
              corr_matrix[,,1] <- corr_celltype1; corr_matrix[,,2] <- corr_celltype2
              for (j in 1:num_celltypes) {
                diag(corr_matrix[,,j]) <- 1 # correlation between the same gene is always one
                cov_tensor[,,j] <- .diag %*% corr_matrix[,,j] %*% .diag  # cov(X) = diag(var(X))^1/2 * corr(X) * diag(var(X))^1/2
              }
              simulated_data <- simulate_bulk_mixture (mu, cov_tensor, p, n=n); Y <- simulated_data$Y
              ##-------------------------
              ##  estimate ratios
              ##-------------------------
              true_theta <- list(p = p, mu = mu, sigma = cov_tensor) %>% enforce_parameter_identifiability()
              overlap <- MixSim::overlap(Pi=p, Mu=t(mu), S=cov_tensor)$BarOmega %>% signif(digits = 3)
              # hellinger <- hellinger_average(p, mu, cov_tensor) %>% round(digits = 4)
              suffix_scenario <- ifelse (.name=="homoscedasctic", "Ho", "He")
              
              estimated_ratios <- suppressWarnings(deconvolute_ratios (mu, Y, scaled=scaled, true_ratios=p, Sigma=cov_tensor,
                                                      deconvolution_functions=deconvolution_functions))
              simulation_metrics_per_config <- tibble::tibble(proportions=name_balance, true_parameters = list(as.list(true_theta)),
                                                              correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
                                                              variance=.name, overlap=overlap, ID = paste0("B", id_scenario, "_", suffix_scenario)) %>%
                dplyr::bind_cols(estimated_ratios)
              
              return(simulation_metrics_per_config)
            }) # end loop scenario variance
          ) 
        }  # end loop correlation second gene
      } # end loop correlation first gene
      simulation_metrics <- simulation_metrics %>%
        dplyr::mutate(entropy=compute_shannon_entropy(p) %>% round(digits = 3))
      dir.create("./simulations/results", showWarnings = F, recursive = TRUE)
      saveRDS(simulation_metrics, file = file.path("./simulations/results", paste0("temp_bivariate_", id_scenario, ".rds")))
      id_scenario <- id_scenario + 1
      return(simulation_metrics)
    }) # end loop ratios
    message("One scenario has been ended.\n\n")
    return(simulations_per_ratios %>% tibble::add_column(centroids=name_distance_centroids))
  }) # end loop mean signatures
  
  return(simulations)
}
