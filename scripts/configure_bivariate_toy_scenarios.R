# Build scenario configuration tibbles for the bivariate (G=2, J=2) toy
# study. Parameter grids live here (not in the package); pass the result to
# run_simulation_benchmark() with a deconvolution_functions list.
#
# Run from the package root after devtools::load_all():
#   source("scripts/configure_bivariate_toy_scenarios.R")
#   cfg <- build_bivariate_scenario_config(corr_sequence = 0)
#   out <- run_simulation_benchmark(cfg, deconvolution_functions, n = 2)

#' Build bivariate generative-model scenario configuration
#'
#' @param proportions Named list of simplex vectors \eqn{\boldsymbol{p}}.
#' @param signature_matrices Named list of mean matrices
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{2\times 2}^{+}}.
#' @param corr_sequence Numeric sequence of within-cell-type correlations.
#' @param diagonal_terms Named list of diagonal variance templates.
#'
#' @return Tibble with `true_theta` list column plus scenario metadata
#'   (`ID`, `overlap`, `entropy`, etc.) for [run_simulation_benchmark()].
build_bivariate_scenario_config <- function(
  proportions = list(
    "balanced" = c(0.5, 0.5),
    "small unbalanced" = c(0.6, 0.4),
    "highly unbalanced" = c(0.05, 0.95)
  ),
  signature_matrices = list(
    "small OVL" = matrix(c(20, 40, 40, 20), nrow = 2)
  ),
  corr_sequence = seq(-0.8, 0.8, 0.2),
  diagonal_terms = list(
    "homoscedastic" = c(1, 1),
    "heteroscedastic" = c(1, 2)
  )
) {
  if (!requireNamespace("MixSim", quietly = TRUE)) {
    stop(
      "build_bivariate_scenario_config() requires MixSim.",
      call. = FALSE
    )
  }

  num_celltypes <- ncol(signature_matrices[[1L]])
  num_genes <- nrow(signature_matrices[[1L]])
  signature_matrices <- purrr::map(
    signature_matrices,
    function(.mean_signature_matrix) {
      dimnames(.mean_signature_matrix) <- list(
        paste0("gene_", seq_len(num_genes)),
        paste0("celltype_", seq_len(num_celltypes))
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

  purrr::pmap(
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
      mu <- signature_matrices[[centroids]]
      p <- proportion_list[[proportion_name]]
      diag_terms <- diagonal_terms[[variance]]

      corr_matrix <- array(
        0,
        dim = c(num_genes, num_genes, num_celltypes),
        dimnames = list(
          paste0("gene_", seq_len(num_genes)),
          paste0("gene_", seq_len(num_genes)),
          paste0("celltype_", seq_len(num_celltypes))
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

      true_theta <- list(p = p, mu = mu, sigma = Sigma)
      overlap_fit <- MixSim::overlap(
        Pi = p,
        Mu = t(mu),
        S = Sigma
      )
      overlap <- signif(overlap_fit$BarOmega, digits = 3)

      tibble::tibble(
        ID = ID,
        correlation_celltype1 = correlation_celltype1,
        correlation_celltype2 = correlation_celltype2,
        overlap = overlap,
        entropy = round(compute_shannon_entropy(p), digits = 3),
        proportions = proportion_name,
        variance = variance,
        centroids = centroids,
        true_theta = list(true_theta)
      )
    }
  ) |>
    dplyr::bind_rows()
}

#' Default deconvolution solvers for the bivariate toy benchmark
#'
#' Omits `lm` / lsfit and CIBERSORT (too few genes for nu-SVR tuning).
#' DeCovarT optimisers use reduced `itmax` for quick smoke runs; increase
#' in production scripts.
#'
#' @param itmax Maximum iterations for gradient-based DeCovarT solvers.
#' @param epsilon Convergence tolerance for DeCovarT solvers.
#'
#' @return Named list suitable for [run_simulation_benchmark()].
bivariate_toy_deconvolution_functions <- function(
  itmax = 10L,
  epsilon = 1e-3
) {
  list(
    "nnls" = list(FUN = deconvolute_ratios_nnls),
    "lsei" = list(FUN = deconvolute_ratios_deconrnaseq),
    "LBFGS" = list(
      FUN = deconvolute_ratios_L_BFGS_B,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "gradient" = list(
      FUN = deconvolute_ratios_gradient_descent,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "Newton-Raphson" = list(
      FUN = deconvolute_ratios_Newton_Raphson,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "Marquardt-Levenberg" = list(
      FUN = deconvolute_ratios_Marquardt_Levenberg,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    ),
    "SA" = list(
      FUN = deconvolute_ratios_simulated_annealing,
      additional_parameters = list(epsilon = epsilon, itmax = itmax)
    )
  )
}
