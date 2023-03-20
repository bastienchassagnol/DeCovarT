#' Simulate a bulk mixture with two cell ratios
#'
#' @param signature_matrix
#' @param cov_matrix
#' @param proportions
#' @param n
#'
#' @return list with X: simulated purified expression profiles (as an array) and Y: the reconstituted bulk mixtures
#' @export


simulate_bulk_mixture <- function(signature_matrix, cov_matrix,
                                  proportions=rep(1/ncol(signature_matrix), ncol(signature_matrix)), n=500) {

  ##################################################################
  ##                        check validity                        ##
  ##################################################################
  if (!all.equal(row.names(signature_matrix), union(dimnames(cov_matrix)[[1]], dimnames(cov_matrix)[[2]])))
    stop("Some of the genes are distinct between expected and covariance expression")
  else if (!all.equal(colnames(signature_matrix), dimnames(cov_matrix)[[3]]))
    stop("Cell types differ between expected and covariance expression")

  ##################################################################
  ##            remove undefined covariance structures            ##
  ##################################################################
  # first remove all cell types associated to a undefinite cov matrix
  valid_celltypes <- purrr::discard(colnames(signature_matrix), ~ cov_matrix %>% magrittr::extract(,,.x) %>% is.na() %>% any())
  # second, remove all cell types associated to a non-positive definite matrix
  valid_celltypes <- purrr::keep(valid_celltypes, ~ cov_matrix %>% magrittr::extract(,,.x) %>% is_positive_definite())
  cov_matrix <- cov_matrix[,,valid_celltypes] # only keep cell_types associated to sound covariance structures
  signature_matrix <- signature_matrix[, valid_celltypes]
  # print(paste("Cell types kept (positive semi-definite covariance matrices) are: ", paste(valid_celltypes, collapse = ", "), "."))

  # initialize Y bulk mixture matrix
  names_genes <- row.names(signature_matrix)
  Y <- matrix(0, nrow = nrow(signature_matrix), ncol = n, dimnames = list(names_genes, paste0("sample_", 1:n)))


  ##### simulation
  X <- array(0, c(nrow(signature_matrix), ncol(signature_matrix),n),
             dimnames = list(names_genes, valid_celltypes, colnames(Y))) # store cell-type specific expression

  for (cell_name in valid_celltypes) { # generate individually each cell distribution
    mean_parameter <- signature_matrix[, cell_name]
    covariance_parameter <- cov_matrix [,,cell_name]
    expression_per_celltype <- MASS::mvrnorm(n = n, mu=mean_parameter, Sigma=covariance_parameter,
                                             tol = 1e-12, empirical = FALSE) # simulate a draw from the covariance matrix
    X[,cell_name, ] <- t(expression_per_celltype)
  }
  for (i in 1:n) {
    Y[,i] <- X[,,i] %*% proportions
  }
  return(list(X=X, Y=Y))
}



# simulate covariance distribution from user-defined estimates
benchmark_deconvolution_algorithms <- function(proportions=list("balanced"=c(0.5, 0.5), "small unbalanced"=c(0.6, 0.4),"highly unbalanced"=c(0.05, 0.95)),
                                               signature_matrix=matrix(c(20, 25, 25, 20), nrow = 2),n=200, score_variable="model_mse", scaled=F,
                                               deconvolution_functions=list("lm" = list(FUN=deconvolute_ratios_abbas,additionnal_parameters=NULL))) {
  ##################################################################
  ##            iterate over covariance structures            ##
  ##################################################################
  num_celltypes <- ncol(signature_matrix); num_genes <- nrow(signature_matrix)
  dimnames(signature_matrix) <- list(paste0("gene_", 1:num_genes), paste0("celltype_", 1:num_celltypes))
  score_variable <- match.arg(score_variable, c("model_mse", "model_rmse", "model_coef_determination", "model_coef_determination_adjusted",
                                                "model_mae", "model_cor", "model_ccc"))

  simulations <- purrr::imap(proportions, function(p, name_scenario) {
    corr_sequence <- seq(-0.8, 0.8, 0.2)
    corr_names <- tidyr::crossing(paste0("correlation_celltype1_", corr_sequence), paste0("_celltype2_", corr_sequence))
    corr_names <- paste0(corr_names[[1]], corr_names[[2]])
    # simulations_storing <- list("homoscedastic"=array(0, c(num_genes,n,length(corr_names)),
    #                                                   dimnames = list(paste0("gene_", 1:num_genes),
    #                                                                   paste0("sample_", 1:n),
    #                                                                   corr_names)),
    #                             "heteroscedastic"=array(0, c(num_genes,n,length(corr_names)),
    #                                                     dimnames = list(paste0("gene_", 1:num_genes),
    #                                                                     paste0("sample_", 1:n),
    #                                                                     corr_names)))
    simulations_dp <- tibble::tibble()
    simulation_metrics <- tibble::tibble()
    for (corr_celltype1 in corr_sequence) {
      print(paste("We are at scenario", name_scenario, "with correlation for celltype 1:", corr_celltype1, "."))
      for (corr_celltype2 in corr_sequence) {
        ##################################################################
        ##                     homoscedastic case                       ##
        ##################################################################
        corr_matrix <- array(0, dim = c(num_genes,num_genes,num_celltypes),
                             dimnames = list(paste0("gene_", 1:num_genes),
                                             paste0("gene_", 1:num_genes),
                                             paste0("celltype_", 1:num_celltypes)))
        cov_matrix <- corr_matrix
        corr_matrix[,,1] <- corr_celltype1; corr_matrix[,,2] <- corr_celltype2
        for (j in 1:num_celltypes) {
          diag(corr_matrix[,,j]) <- 1 # correlation between the same gene is always one
          cov_matrix[,,j] <- diag(c(1, 1)) %*% corr_matrix[,,j] %*% diag(c(1, 1))  # cov(X) = diag(var(X))^1/2 * corr(X) * diag(var(X))^1/2
        }
        overlap_homoscedastic <- MixSim::overlapGOM(Pi = p, Mu = signature_matrix, S = cov_matrix) %>% round(digits = 3)
        hellinger_homoscedastic <- hellinger_average(p, signature_matrix, cov_matrix) %>% round(digits = 4)


        simulated_data <- simulate_bulk_mixture (signature_matrix, cov_matrix, p, n=n); Y <- simulated_data$Y
        title_dp <- paste0("Correlation celltype 1:", corr_celltype1, "celltype 2:", corr_celltype2)
        dp <- plot_multivariate_density(simulated_data$X) +
          ggtitle(paste(title_dp, "in homo scenario"), subtitle = paste("overlap:", overlap_homoscedastic,
                                                                        "hellinger:", hellinger_homoscedastic))


        estimated_ratios <- deconvolute_ratios (signature_matrix, Y, scaled=F, true_ratios=p, Sigma=cov_matrix,
                                                deconvolution_functions=deconvolution_functions)
        # simulations_storing[["homoscedastic"]][,,title_dp] <- Y
        simulation_metrics <- simulation_metrics %>% dplyr::bind_rows(
          tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
                         variance="homoscedastic", overlap=overlap_homoscedastic, hellinger=hellinger_homoscedastic) %>%
            dplyr::bind_cols(estimated_ratios))
        simulations_dp <- simulations_dp %>% dplyr::bind_rows(
          tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
                         variance="homoscedastic", overlap=overlap_homoscedastic, hellinger=hellinger_homoscedastic, density_plot=list(dp)))
        ##################################################################
        ##                     heteroscedastic case                     ##
        ##################################################################

        for (j in 1:num_celltypes) {
          cov_matrix[,,j] <- diag(c(1, 4)) %*% corr_matrix[,,j] %*% diag(c(1, 4))  # some genes have variance 1, others variance 2
        }
        overlap_heteroscedastic <- MixSim::overlapGOM(Pi = p, Mu = signature_matrix, S = cov_matrix) %>% round(digits = 3)
        hellinger_heteroscedastic <- hellinger_average(p, signature_matrix, cov_matrix) %>% round(digits = 4)

        simulated_data <- simulate_bulk_mixture (signature_matrix, cov_matrix, p, n=n); Y <- simulated_data$Y
        dp <-  plot_multivariate_density(simulated_data$X) +
          ggtitle(paste(title_dp, "in hetero scenario"), subtitle = paste("overlap:", overlap_heteroscedastic,
                                                                          "hellinger:", hellinger_heteroscedastic))

        estimated_ratios <- deconvolute_ratios (signature_matrix, Y, scaled=F, true_ratios=p, Sigma=cov_matrix,
                                                deconvolution_functions=deconvolution_functions)
        # simulations_storing[["heteroscedastic"]][,,title_dp] <- Y
        simulation_metrics <- simulation_metrics %>% dplyr::bind_rows(
          tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1,
                         overlap=overlap_heteroscedastic, hellinger=hellinger_heteroscedastic,
                         correlation_celltype2=corr_celltype2,variance="heteroscedastic") %>%
            dplyr::bind_cols(estimated_ratios))
        simulations_dp <- simulations_dp %>% dplyr::bind_rows(
          tibble::tibble(proportions=name_scenario, correlation_celltype1=corr_celltype1, correlation_celltype2=corr_celltype2,
                         variance="heteroscedastic", overlap=overlap_heteroscedastic, hellinger=hellinger_heteroscedastic, density_plot=list(dp)))

      }  # end loop correlation second gene
    } # end loop correlation first gene
    entropy <- compute_shannon_entropy(p) %>% round(digits = 3); simulation_metrics <- simulation_metrics %>% dplyr::mutate(entropy=entropy)

    ##################################################################
    ##                       Complex Heatmaps                       ##
    ##################################################################
    CH_homo <- purrr::imap(split(simulation_metrics %>% filter(variance=="homoscedastic"), simulation_metrics$deconvolution_name), function(.x, .y) {
      cor_matrix_per_algo <- .x %>% select (all_of(c("correlation_celltype1", "correlation_celltype2", score_variable))) %>%
        tidyr::pivot_wider(names_from = c(correlation_celltype2), values_from = .data[[score_variable]],values_fn = mean) %>%
        tibble::column_to_rownames("correlation_celltype1") %>% as.matrix()

      complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(cor_matrix_per_algo, name=gsub("model_", "", score_variable),
                                                          heatmap_legend_param = list(title = gsub("model_", "", score_variable)), row_title="Correlation cell type 1",
                                                          cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix_per_algo),
                                                          row_title_gp = grid::gpar(fontsize = 10), column_names_rot = 45,
                                                          cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
                                                          column_labels = colnames(cor_matrix_per_algo),
                                                          width = unit(8, "cm"), height = unit(8, "cm"),
                                                          column_title_gp = grid::gpar(fontsize = 10), column_title="Correlation cell type 2") %>%
        ComplexHeatmap::draw(padding=unit(c(0, 0, 0, 0), "cm"),
                             column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>%
        grid::grid.grabExpr()
      return(complex_heatmap_per_algo)
    })
    CH_homo <- gridExtra::arrangeGrob(grobs = CH_homo, ncol = min(length(deconvolution_functions), 2),  padding = unit(0.1, "line"),
                                      top=ggpubr::text_grob(paste("Entropy is:", entropy, "\n homo scenario"), size = 18, face = "bold"))

    CH_hetero <- purrr::imap(split(simulation_metrics %>% filter(variance=="heteroscedastic"), simulation_metrics$deconvolution_name), function(.x, .y) {
      cor_matrix_per_algo <- .x %>% select (all_of(c("correlation_celltype1", "correlation_celltype2", score_variable))) %>%
        tidyr::pivot_wider(names_from = c(correlation_celltype2), values_from = .data[[score_variable]],values_fn = mean) %>%
        tibble::column_to_rownames("correlation_celltype1") %>% as.matrix()

      complex_heatmap_per_algo <- ComplexHeatmap::Heatmap(cor_matrix_per_algo, name=gsub("model_", "", score_variable),
                                                          heatmap_legend_param = list(title = gsub("model_", "", score_variable)), row_title="Correlation cell type 1",
                                                          cluster_rows = F, row_names_gp = grid::gpar(fontsize = 8), row_labels = colnames(cor_matrix_per_algo),
                                                          row_title_gp = grid::gpar(fontsize = 10), column_names_rot = 45,
                                                          cluster_columns = F, column_names_gp = grid::gpar(fontsize = 8),
                                                          column_labels = colnames(cor_matrix_per_algo),
                                                          width = unit(8, "cm"), height = unit(8, "cm"),
                                                          column_title_gp = grid::gpar(fontsize = 10), column_title="Correlation cell type 2") %>%
        ComplexHeatmap::draw(padding=unit(c(0, 0, 0, 0), "cm"),
                             column_title= .y, column_title_gp = grid::gpar(fontsize = 12, fontface="bold")) %>%
        grid::grid.grabExpr()
      return(complex_heatmap_per_algo)
    })
    CH_hetero <- gridExtra::arrangeGrob(grobs = CH_hetero, ncol = min(length(deconvolution_functions), 2),  padding = unit(0.1, "line"),
                                        top=ggpubr::text_grob(paste("Entropy is:", entropy, "\n hetero scenario"), size = 18, face = "bold"))


    #################################################################
    ##                        mse boxplots                        ##
    #################################################################

    mse_boxplots <- plot_boxplots_parameters(simulation_metrics, score_variable) +
      ggtitle(paste("Entropy is:", entropy))

    metric_plots <- gridExtra::arrangeGrob(grobs=list(CH_homo, CH_hetero, mse_boxplots),
                                           layout_matrix = matrix(c(1, 2, 3, 3), nrow = 2, byrow = T), heights = c(1.5, 2))

    return(list("simulation_metrics"=simulation_metrics,
                "metric_plots"=metric_plots, "simulations_dp"=simulations_dp))
  })

  return(simulations)
}























