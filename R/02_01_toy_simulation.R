#' Simulate bulk mixtures from a Gaussian convolution
#'
#' @description
#' Draws purified profiles
#' \eqn{\boldsymbol{x}_j\sim\mathcal{N}_{G}(\boldsymbol{\mu}_{\cdot j},
#' \boldsymbol{\Sigma}_j)} independently for each cell type and forms bulk
#' observations by the linear mixture
#' \deqn{
#'   \boldsymbol{y}_{\cdot i}
#'   =\boldsymbol{\mu}^{(i)}\boldsymbol{p}
#'   =\sum_{j=1}^{J}p_j\,\boldsymbol{x}_j^{(i)},
#' }
#' matching the article's conditional model
#' \eqn{\boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})\sim
#' \mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \boldsymbol{\Sigma}(\boldsymbol{p}))} with
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^{2}\boldsymbol{\Sigma}_j}.
#'
#' @param signature_matrix Mean matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}.
#' @param Sigma Array of covariances
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}}.
#' @param p Proportion vector \eqn{\boldsymbol{p}\in\Delta^{J-1}}
#'   (default: uniform).
#' @param n Number of bulk samples \eqn{n}.
#'
#' @return A list with:
#' * `mean_signature_matrix`: array
#'   \eqn{(x_{gji})\in\mathcal{M}_{G\times J\times N}} of simulated purified
#'   profiles;
#' * `Y`: matrix \eqn{\boldsymbol{Y}\in\mathcal{M}_{G\times N}} whose columns
#'   are bulk vectors \eqn{\boldsymbol{y}_{\cdot i}}.
#'
#' @export
#' @seealso [deconvolute_ratios()], [benchmark_bivariate_gaussian_convolutions()]
simulate_bulk_mixture <- function(
  signature_matrix,
  Sigma,
  p = rep(1 / ncol(signature_matrix), ncol(signature_matrix)),
  n = 500
) {
  ##################################################################
  ##                        check validity                        ##
  ##################################################################
  if (
    !all.equal(
      row.names(signature_matrix),
      union(dimnames(Sigma)[[1]], dimnames(Sigma)[[2]])
    )
  ) {
    stop(
      "Some of the genes are distinct between expected and covariance expression"
    )
  } else if (!all.equal(colnames(signature_matrix), dimnames(Sigma)[[3]])) {
    stop("Cell types differ between expected and covariance expression")
  }

  ##### simulation
  valid_celltypes <- colnames(signature_matrix)
  names_genes <- row.names(signature_matrix)
  Y <- matrix(
    0,
    nrow = nrow(signature_matrix),
    ncol = n,
    dimnames = list(names_genes, paste0("sample_", 1:n))
  )
  mean_signature_matrix <- array(
    0,
    c(nrow(signature_matrix), ncol(signature_matrix), n),
    dimnames = list(names_genes, valid_celltypes, colnames(Y))
  ) # store cell-type specific expression

  for (cell_name in valid_celltypes) {
    # generate individually each cell distribution
    mean_parameter <- signature_matrix[, cell_name]
    covariance_parameter <- Sigma[,, cell_name]
    expression_per_celltype <- MASS::mvrnorm(
      n = n,
      mu = mean_parameter,
      Sigma = covariance_parameter,
      tol = 1e-12,
      empirical = FALSE
    ) # simulate a draw from the covariance matrix
    mean_signature_matrix[, cell_name, ] <- t(expression_per_celltype)
  }

  Y <- tensor::tensor(p, B = mean_signature_matrix, alongA = 1, alongB = 2) # tensor product does directly the computation mean_signature_matrix %*% p
  return(list(mean_signature_matrix = mean_signature_matrix, Y = Y))
}


#' Benchmark bivariate Gaussian convolutions
#'
#' @description
#' Wrapper reproducing the bivariate (\eqn{G=2}, \eqn{J=2}) toy study of the
#' article: for each scenario it builds
#' \eqn{\boldsymbol{\mu}} and
#' \eqn{(\boldsymbol{\Sigma}_j)_{j}}, simulates
#' \eqn{\boldsymbol{Y}} via [simulate_bulk_mixture()], and deconvolves with the
#' supplied algorithms. Performance is summarised against entropy of
#' \eqn{\boldsymbol{p}} and overlap of the Gaussian mixture.
#'
#' @details
#' Designed for two cell types and two genes. Larger \eqn{(G,J)} with only
#' bivariate observations is prone to non-identifiability.
#'
#' Scenarios are enumerated with [tidyr::expand_grid()] and tagged with a
#' unique `ID` via [dplyr::row_number()], so no side-effect mutation is needed
#' while looping over the design.
#'
#' @param proportions List of simplex vectors \eqn{\boldsymbol{p}}.
#' @param signature_matrices List of mean matrices
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{2\times 2}^{+}}.
#' @param corr_sequence,diagonal_terms Correlation sequence and diagonal
#'   variance templates used to assemble
#'   \eqn{\boldsymbol{\Sigma}_j=\mathrm{D}_{j}^{1/2}\mathbf{R}_j\mathrm{D}_{j}^{1/2}}.
#' @param deconvolution_functions Named list of deconvolution callables
#'   (each with `FUN` and optional `additional_parameters`).
#' @param n Number of bulk replicates \eqn{N} per scenario.
#' @param scaled Logical; whether to log-scale inputs before estimation.
#' @param cores Parallel workers for [deconvolute_ratios()].
#'
#' @return A list with `config` (design + entropy/overlap) and `simulations`
#'   (estimation tibble).
#'
#' @importFrom rlang .data
#' @export
#' @seealso [simulate_bulk_mixture()], [deconvolute_ratios()]
benchmark_bivariate_gaussian_convolutions <- function(
  proportions = list(
    "balanced" = c(0.5, 0.5),
    "small unbalanced" = c(0.6, 0.4),
    "highly unbalanced" = c(0.05, 0.95)
  ),
  signature_matrices = list("small OVL" = matrix(c(20, 40, 40, 20), nrow = 2)),
  corr_sequence = seq(-0.8, 0.8, 0.2),
  diagonal_terms = list("homoscedastic" = c(1, 1), "heteroscedastic" = c(1, 2)),
  deconvolution_functions = list(
    "lsfit" = list(FUN = deconvolute_ratios_lsfit, additional_parameters = NULL)
  ),
  n = 200,
  scaled = FALSE,
  cores = ifelse(
    .Platform$OS.type == "unix",
    getOption("mc.cores", parallel::detectCores()),
    1
  )
) {
  num_celltypes <- ncol(signature_matrices[[1]])
  num_genes <- nrow(signature_matrices[[1]])
  signature_matrices <- purrr::map(
    signature_matrices,
    function(.mean_signature_matrix) {
      dimnames(.mean_signature_matrix) <- list(
        paste0("gene_", 1:num_genes),
        paste0("celltype_", 1:num_celltypes)
      )
      .mean_signature_matrix
    }
  )

  proportion_list <- proportions
  design <- tidyr::expand_grid(
    centroids = names(signature_matrices),
    proportion_name = names(proportion_list),
    correlation_celltype1 = corr_sequence,
    correlation_celltype2 = corr_sequence,
    variance = names(diagonal_terms)
  ) |>
    dplyr::mutate(
      scenario_idx = dplyr::row_number(),
      ID = paste0(
        "B",
        .data$scenario_idx,
        "_",
        ifelse(.data$variance == "homoscedastic", "Ho", "He")
      )
    )

  scenario_results <- purrr::pmap(
    design,
    function(
      centroids,
      proportion_name,
      correlation_celltype1,
      correlation_celltype2,
      variance,
      scenario_idx,
      ID
    ) {
      message(paste0(
        "Scenario ", ID, ": ", proportion_name,
        ", corr=(", correlation_celltype1, ", ",
        correlation_celltype2, "), centroids=", centroids,
        ", variance=", variance, "."
      ))

      mu <- signature_matrices[[centroids]]
      p <- proportion_list[[proportion_name]]
      diag_terms <- diagonal_terms[[variance]]

      corr_matrix <- array(
        0,
        dim = c(num_genes, num_genes, num_celltypes),
        dimnames = list(
          paste0("gene_", 1:num_genes),
          paste0("gene_", 1:num_genes),
          paste0("celltype_", 1:num_celltypes)
        )
      )
      Sigma <- corr_matrix
      corr_matrix[,, 1] <- correlation_celltype1
      corr_matrix[,, 2] <- correlation_celltype2
      for (j in seq_len(num_celltypes)) {
        diag(corr_matrix[,, j]) <- 1
        Sigma[,, j] <- sqrt(diag(diag_terms)) %*%
          corr_matrix[,, j] %*%
          sqrt(diag(diag_terms))
      }

      simulated_data <- simulate_bulk_mixture(
        signature_matrix = mu,
        Sigma = Sigma,
        p = p,
        n = n
      )
      true_theta <- list(p = p, mu = mu, sigma = Sigma) |>
        enforce_parameter_identifiability()
      overlap <- MixSim::overlap(
        Pi = p,
        Mu = t(mu),
        S = Sigma
      )$BarOmega |>
        signif(digits = 3)

      estimated_ratios <- suppressWarnings(deconvolute_ratios(
        signature_matrix = mu,
        bulk_expression = simulated_data$Y,
        true_ratios = p,
        Sigma = Sigma,
        deconvolution_functions = deconvolution_functions,
        scaled = scaled,
        cores = cores
      ))

      simulations <- tibble::tibble(
        ID = ID,
        correlation_celltype1 = correlation_celltype1,
        correlation_celltype2 = correlation_celltype2
      ) |>
        dplyr::bind_cols(estimated_ratios)

      config <- tibble::tibble(
        ID = ID,
        overlap = overlap,
        entropy = round(compute_shannon_entropy(p), digits = 3),
        proportions = proportion_name,
        variance = variance,
        centroids = centroids,
        true_parameters = list(as.list(true_theta)),
        nobservations = n
      )

      list(simulations = simulations, config = config)
    }
  )

  list(
    simulations = dplyr::bind_rows(purrr::map(scenario_results, "simulations")),
    config = dplyr::bind_rows(purrr::map(scenario_results, "config"))
  )
}
