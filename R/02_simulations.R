#' Simulate bulk mixtures
#'
#' @description Using user defined parameters \eqn{p}, \eqn{\boldsymbol{X}} and
#' \eqn{\boldsymbol{cov}}, simulate virtual bulk mixtures, using standard linear
#' constraint of deconvolution algorithms, namely:
#' \deqn{\boldsymbol{\hat{y}}=\boldsymbol{X} \times \hat{\boldsymbol{p}}}
#'
#' @inheritParams deconvolute_ratios
#' @param p The ratios to estimate (by default, equi-balanced cellular
#' proportions are simulated)
#' @param n The integer number of samples to generate
#'
#' @import ggplot2
#' @return list with two items:
#' * \eqn{\boldsymbol{X}:\left(x_{g,j,i}\right)_{1\le g\le G, 1\le j \le J, 1\le i \le N}}:
#' the simulated purified expression profiles, stored in
#' a three-dimensional array, each array storing the simulated expression of the
#' \eqn{J} cell populations for a given individual \eqn{i}
#' * \eqn{\boldsymbol{Y} \in \mathcal{M}_{G \times N}^+}: a two-dimensional matrix storing the virtually
#' simulated bulk mixtures. Each column \eqn{\boldsymbol{y_{.i}}} stores for individual \eqn{i}
#' the expression value of the \eqn{G} genes, while attribute row.names is used
#' to uniquely identify each of them.
#' @export



simulate_bulk_mixture <- function(signature_matrix, Sigma,
                                  p = rep(1 / ncol(signature_matrix), ncol(signature_matrix)), n = 500) {
  ##################################################################
  ##                        check validity                        ##
  ##################################################################
  if (!all.equal(row.names(signature_matrix), union(dimnames(Sigma)[[1]], dimnames(Sigma)[[2]]))) {
    stop("Some of the genes are distinct between expected and covariance expression")
  } else if (!all.equal(colnames(signature_matrix), dimnames(Sigma)[[3]])) {
    stop("Cell types differ between expected and covariance expression")
  }

  ##################################################################
  ##            remove undefined covariance structures            ##
  ##################################################################
  # first remove all cell types associated to a undefinite cov matrix
  valid_celltypes <- purrr::discard(colnames(signature_matrix), ~ Sigma %>%
    magrittr::extract(, , .x) %>%
    is.na() %>%
    any())
  # second, remove all cell types associated to a non-positive definite matrix
  valid_celltypes <- purrr::keep(valid_celltypes, ~ Sigma %>%
    magrittr::extract(, , .x) %>%
    is_positive_definite())
  # only keep cell_types associated to sound covariance structures
  Sigma <- Sigma[, , valid_celltypes]
  signature_matrix <- signature_matrix[, valid_celltypes]

  ##### simulation
  valid_celltypes <- colnames(signature_matrix)
  names_genes <- row.names(signature_matrix)
  Y <- matrix(0, nrow = nrow(signature_matrix), ncol = n, dimnames = list(names_genes, paste0("sample_", 1:n)))
  X <- array(0, c(nrow(signature_matrix), ncol(signature_matrix), n),
    dimnames = list(names_genes, valid_celltypes, colnames(Y))
  ) # store cell-type specific expression

  for (cell_name in valid_celltypes) { # generate individually each cell distribution
    mean_parameter <- signature_matrix[, cell_name]
    covariance_parameter <- Sigma[, , cell_name]
    expression_per_celltype <- MASS::mvrnorm(
      n = n, mu = mean_parameter, Sigma = covariance_parameter,
      tol = 1e-12, empirical = FALSE
    ) # simulate a draw from the covariance matrix
    X[, cell_name, ] <- t(expression_per_celltype)
  }

  # Y_test <- matrix(0, nrow = nrow(signature_matrix), ncol = n, dimnames = list(names_genes, paste0("sample_", 1:n)))
  # for (i in 1:n) {
  #   Y_test[,i] <- X[,,i] %*% proportions
  # }

  Y <- tensor::tensor(p, B = X, alongA = 1, alongB = 2) # tensor product does directly the computation X %*% p
  return(list(X = X, Y = Y))
}




#' Benchmark two-component cell populations
#'
#' @description It is a highly specific wrapper function, used to reproduce the results
#' for the bivariate toy example provided in the paper associated to this package.
#' Precisely, it has been designed to predict the performance of a variety of deconvolution
#' algorithms, in the specific framework where two cell populations are characterised
#' by two genes.
#'
#' @details Purpose of this toy example was to evaluate the quality of the benchmark
#' with respect to the level of entropy (the global population disequilibrium)
#' and the overlap, in a two dimensional convolution of Gaussian mixture distributions.
#'
#'
#' @note In a near future, we would provide an additional and more flexible
#' benchmark, however, till now, we strongly do not recommend to benchmark
#' virtual mixtures with more than two populations, since, with only two genes,
#' you are likely to come across strong identifiability issues.
#'
#' @inheritParams deconvolute_ratios
#' @param proportions A list of numerical vectors, assumed each to respect the unit
#' simplex constraint:
#' \deqn{\begin{cases} \sum_{j=1}^J p_{j}=1\\ \forall j \in \widetilde{J} \quad p_j\ge 0 \end{cases}}
#' @param signature_matrices A list of bivariate matrices, in \eqn{\mathcal{M}_{G \times J}^+}
#'  storing the purified expression profiles of the \eqn{J} cell populations of the sample.
#'  This parameter is especially useful to assert the performance of the deconvolution
#'  algorithms, as a function of the proximity of the respective centroids (mean
#'  expression profiles) of the evaluated cell populations.
#' @param diagonal_terms,corr_sequence Together, these parameters control the
#' characteristics of the array of covariance matrices used to model the transcriptomic
#' network structure of each cell population:
#' * A numeric vector, strictly included between -1 and 1, is expected for
#' param `corr_sequence`. Each paired arrangement of correlation for both cell
#' populations will be used to describe the pairwise correlation level between
#' gene 1 and gene 2.
#' * A list of numeric vectors, used to parametrise the diagonal terms of
#' the covariance matrix, is expected for param `diagonal_terms`. Explicit
#' names can be given as well to properly identify a scenario.  By default,
#' we compare a standard \emph{homoscedastic} scenario, in which we assume a constant diagonal term
#' of 1 for both genes, with a \emph{heteroscedastic} scenario, in which we assign
#' a variance of 2 for gene 2, for both cell populations. Finally, the covariance
#' matrix is computed and normalised back to integrate both constraints, using
#' following matricial relationship, linking the value of a correlation matrix with
#' its respective covariance matrix:
#' \deqn{\text{Cov}(\boldsymbol{X}) = \text{diag}\left[(\text{Var}(\boldsymbol{X})\right]^{1/2} \times \text{Cor}(\boldsymbol{X}) \times \text{diag}\left[(\text{Var}(\boldsymbol{X})\right]^{1/2}}
#' @param n The number of samples to generate (an integer, but upon request, we may
#' update the function to accept a numerical vector instead)
#'
#' @return a two-entry list, with:
#' * `config`: the complete parameter configuration used to generate the benchmark
#' including the list of parameters, as well as some general metrics, such as the
#' overlap and the entropy of the convolution of GMMs used
#' * `simulations`: the results of the simulations themselves,
#' mostly similar to the output of [deconvolute_ratios]
#'
#' @seealso [deconvolute_ratios()]
#' @export


benchmark_deconvolution_algorithms_two_genes <- function(proportions = list("balanced" = c(0.5, 0.5), "small unbalanced" = c(0.6, 0.4), "highly unbalanced" = c(0.05, 0.95)),
                                                         signature_matrices = list("small OVL" = matrix(c(20, 40, 40, 20), nrow = 2)), corr_sequence = seq(-0.8, 0.8, 0.2),
                                                         diagonal_terms = list("homoscedastic" = c(1, 1), "heteroscedastic" = c(1, 2)),
                                                         deconvolution_functions = list("lm" = list(FUN = deconvolute_ratios_abbas, additional_parameters = NULL)),
                                                         n = 200, scaled = FALSE,
                                                         cores = ifelse(.Platform$OS.type == "unix", getOption("mc.cores", parallel::detectCores()), 1)) {
  ##################################################################
  ##            iterate over covariance structures            ##
  ##################################################################
  num_celltypes <- ncol(signature_matrices[[1]])
  num_genes <- nrow(signature_matrices[[1]])
  signature_matrices <- purrr::map(signature_matrices, function(.x) {
    dimnames(.x) <- list(paste0("gene_", 1:num_genes), paste0("celltype_", 1:num_celltypes))
    return(.x)
  })
  id_scenario <- 1
  id_tibble <- tibble::tibble() # initiate ID scenario
  simulations <- purrr::imap_dfr(signature_matrices, function(mu, name_distance_centroids) {
    simulations_per_ratios <- purrr::imap_dfr(proportions, function(p, name_balance) {
      simulation_metrics <- tibble::tibble()
      for (corr_celltype1 in corr_sequence) {
        message(paste0(
          "We are at equilibrium scenario: ", name_balance, " with correlation for celltype 1: ", corr_celltype1,
          " and distance to centroids: ", name_distance_centroids, "."
        ))
        for (corr_celltype2 in corr_sequence) {
          simulation_metrics <- simulation_metrics %>% dplyr::bind_rows(
            purrr::imap_dfr(diagonal_terms, function(.diag, .name) {
              ## ------------------------------
              ##  generate covariance matrix
              ## ------------------------------
              corr_matrix <- array(0,
                dim = c(num_genes, num_genes, num_celltypes),
                dimnames = list(
                  paste0("gene_", 1:num_genes),
                  paste0("gene_", 1:num_genes),
                  paste0("celltype_", 1:num_celltypes)
                )
              )
              Sigma <- corr_matrix
              corr_matrix[, , 1] <- corr_celltype1
              corr_matrix[, , 2] <- corr_celltype2
              for (j in 1:num_celltypes) {
                diag(corr_matrix[, , j]) <- 1 # correlation between the same gene is always one
                Sigma[, , j] <- sqrt(diag(.diag)) %*% corr_matrix[, , j] %*% sqrt(diag(.diag)) # cov(X) = diag(var(X))^1/2 * corr(X) * diag(var(X))^1/2
              }
              simulated_data <- simulate_bulk_mixture(signature_matrix = mu, Sigma = Sigma, p = p, n = n)
              Y <- simulated_data$Y
              ## -------------------------
              ##  estimate ratios
              ## -------------------------
              true_theta <- list(p = p, mu = mu, sigma = Sigma) %>% enforce_parameter_identifiability()
              overlap <- MixSim::overlap(Pi = p, Mu = t(mu), S = Sigma)$BarOmega %>% signif(digits = 3)
              # hellinger <- hellinger_average(p, mu, Sigma) %>% round(digits = 4)

              estimated_ratios <- suppressWarnings(deconvolute_ratios(
                signature_matrix = mu, bulk_expression = Y,
                true_ratios = p, Sigma = Sigma,
                deconvolution_functions = deconvolution_functions,
                scaled = scaled, cores = cores
              ))
              name_id <- paste0("B", id_scenario, "_", ifelse(.name == "homoscedastic", "Ho", "He"))
              simulation_metrics_per_config <- tibble::tibble(
                ID = name_id, correlation_celltype1 = corr_celltype1,
                correlation_celltype2 = corr_celltype2
              ) %>%
                dplyr::bind_cols(estimated_ratios)
              id_tibble_temp <- tibble::tibble(
                ID = name_id,
                overlap = overlap, entropy = compute_shannon_entropy(p) %>% round(digits = 3),
                proportions = name_balance, variance = .name, centroids = name_distance_centroids,
                true_parameters = list(as.list(true_theta)), nobservations = n
              )
              id_tibble <<- id_tibble %>% dplyr::bind_rows(id_tibble_temp)

              return(simulation_metrics_per_config)
            }) # end loop scenario variance
          )
        } # end loop correlation second gene
      } # end loop correlation first gene
      dir.create("./simulations/results", showWarnings = F, recursive = TRUE)
      saveRDS(simulation_metrics, file = file.path("./simulations/results", paste0("temp_bivariate_", id_scenario, ".rds")))

      id_scenario <<- id_scenario + 1
      message("One scenario has been ended.\n\n")
      return(simulation_metrics)
    }) # end loop ratios
    return(simulations_per_ratios)
  }) # end loop mean signatures

  return(list("simulations" = simulations, "config" = id_tibble))
}
