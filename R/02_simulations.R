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
  # for (i in 1:n) {
  #   Y[,i] <- X[,,i] %*% proportions
  # }
  Y <- tensor::tensor(proportions, B = X, alongA = 1, alongB = 2) # tensor product does directly the computation X %*% p
  return(list(X=X, Y=Y))
}




#' Benchmark two-component cell populations
#'
#' @param proportions
#' @param signature_matrix
#' @param n
#' @param score_variable
#' @param scaled
#' @param deconvolution_functions
#'
#' @return
#' @export


benchmark_deconvolution_algorithms <- function(proportions=list("balanced"=c(0.5, 0.5), "small unbalanced"=c(0.6, 0.4),"highly unbalanced"=c(0.05, 0.95)),
                                               signature_matrix=matrix(c(20, 25, 25, 20), nrow = 2), corr_sequence = seq(-0.8, 0.8, 0.2),
                                               deconvolution_functions=list("lm" = list(FUN=deconvolute_ratios_abbas,additionnal_parameters=NULL)),
                                               n=200, scaled=FALSE) {
  ##################################################################
  ##            iterate over covariance structures            ##
  ##################################################################
  num_celltypes <- ncol(signature_matrix); num_genes <- nrow(signature_matrix)
  dimnames(signature_matrix) <- list(paste0("gene_", 1:num_genes), paste0("celltype_", 1:num_celltypes))

  simulations <- purrr::imap_dfr(proportions, function(p, name_scenario) {
    corr_names <- tidyr::crossing(paste0("correlation_celltype1_", corr_sequence), paste0("_celltype2_", corr_sequence))
    corr_names <- paste0(corr_names[[1]], corr_names[[2]])
    simulations_dp <- tibble::tibble();     simulation_metrics <- tibble::tibble()
    for (corr_celltype1 in corr_sequence) {
      print(paste("We are at scenario", name_scenario, "with correlation for celltype 1:", corr_celltype1, "."))
      for (corr_celltype2 in corr_sequence) {
        ##################################################################
        ##                     homoscedastic case                       ##
        ##################################################################
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
          cov_tensor[,,j] <- diag(c(1, 1)) %*% corr_matrix[,,j] %*% diag(c(1, 1))  # cov(X) = diag(var(X))^1/2 * corr(X) * diag(var(X))^1/2
        }
        simulated_data <- simulate_bulk_mixture (signature_matrix, cov_tensor, p, n=n); Y <- simulated_data$Y
        ##-------------------------
        ##  distribution features
        ##-------------------------
        overlap_homoscedastic <- compute_average_overlap (list(p=p, mu=signature_matrix,sigma=cov_tensor)) %>% round(digits = 3)
        hellinger_homoscedastic <- hellinger_average(p, signature_matrix, cov_tensor) %>% round(digits = 4)


        estimated_ratios <- deconvolute_ratios (signature_matrix, Y, scaled=scaled, true_ratios=p, Sigma=cov_tensor,
                                                deconvolution_functions=deconvolution_functions)
        simulation_metrics <- simulation_metrics %>% dplyr::bind_rows(
          tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
                         variance="homoscedastic", overlap=overlap_homoscedastic, hellinger=hellinger_homoscedastic) %>%
            dplyr::bind_cols(estimated_ratios))
        ##################################################################
        ##                     heteroscedastic scenario         ##
        ##################################################################

        for (j in 1:num_celltypes) {
          cov_tensor[,,j] <- diag(c(1, sqrt(2))) %*% corr_matrix[,,j] %*% diag(c(1, sqrt(2)))  # some genes have var 1, others var 2
        }
        overlap_heteroscedastic <- MixSim::overlapGOM(Pi = p, Mu = signature_matrix, S = cov_tensor) %>% round(digits = 3)
        hellinger_heteroscedastic <- hellinger_average(p, signature_matrix, cov_tensor) %>% round(digits = 4)

        simulated_data <- simulate_bulk_mixture (signature_matrix, cov_tensor, p, n=n); Y <- simulated_data$Y
        estimated_ratios <- deconvolute_ratios (signature_matrix, Y, scaled=F, true_ratios=p, Sigma=cov_tensor,
                                                deconvolution_functions=deconvolution_functions)
        simulation_metrics <- simulation_metrics %>% dplyr::bind_rows(
          tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1,
                         overlap=overlap_heteroscedastic, hellinger=hellinger_heteroscedastic,
                         correlation_celltype2=corr_celltype2,variance="heteroscedastic") %>%
            dplyr::bind_cols(estimated_ratios))

      }  # end loop correlation second gene
    } # end loop correlation first gene
    entropy <- compute_shannon_entropy(p) %>% round(digits = 3)
    simulation_metrics <- simulation_metrics %>% dplyr::mutate(entropy=entropy)
    return(simulation_metrics)
  })

  return(simulations)
}
