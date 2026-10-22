# Single-cell reference moments.
# The network generator shared by both models is
# R/01_02_generate_synthetic_networks.R.
# This file only rescales those moments onto a cell-size convention.

#' Cell-size conventions for a single-cell reference
#'
#' @description
#' Turns nonnegative gene-by-cell-type moments into the reference used
#' by the single-cell mixture. Column \eqn{j} of the supplied mean is
#' normalised to a composition \eqn{\boldsymbol{\theta}_{\cdot j}}. Two
#' conventions are available; the default is the raw-covariance,
#' cell-fraction model of the article.
#'
#' * `cellular_raw` (default). Mean
#'   \eqn{\tilde{\boldsymbol{\mu}}_{\cdot j}=S_j\boldsymbol{\theta}_{\cdot j}}
#'   with the covariance left as supplied. Keeping \eqn{S_j} in the
#'   mean makes the fitted composition a cell fraction, not a
#'   transcript fraction.
#' * `transcript_raw`. Library-normalised mean
#'   \eqn{\boldsymbol{\theta}_{\cdot j}} and per-transcript covariance
#'   \eqn{\boldsymbol{\Sigma}_j/S_j}. The fitted composition is the
#'   transcript fraction \eqn{\boldsymbol{q}}. The physical generator
#'   is still the cellular model: `generative_mean` and
#'   `generative_sigma` are the moments [simulate_bulk_from_cells()]
#'   should draw from.
#'
#' @param moments List with `mean_profiles` (\eqn{G\times J}) and
#'   `covariance_matrices` (\eqn{G\times G\times J}), as returned by
#'   [simulate_hierarchical_grn_moments()] with `nonnegative = TRUE`.
#' @param cell_size Positive cell sizes \eqn{S_j}. A scalar is repeated
#'   once per cell type. For two types the scenarios use \eqn{(1, 1)}
#'   and \eqn{(1, 5)}; \eqn{(1, 2, 5)} is the same idea for three types.
#' @param reference `cellular_raw` or `transcript_raw`.
#'
#' @return A list with the matrices passed to the estimator
#'   (`mean_signature_matrix`, `Sigma`), the cellular generator
#'   (`generative_mean`, `generative_sigma`), `cell_size`,
#'   `composition` (\eqn{\boldsymbol{\theta}}), `reference` and
#'   `estimand` (`"cellular"` or `"transcript"`).
#'
#' @examples
#' mu <- cbind(c(2, 2, 0, 0), c(0, 0, 1, 1))
#' sigma <- array(rep(diag(4), 2), dim = c(4, 4, 2))
#' moments <- list(mean_profiles = mu, covariance_matrices = sigma)
#' ref <- rescale_single_cell_moments(moments, c(1, 5))
#' colSums(ref$generative_mean)
#' @export
#' @seealso [simulate_hierarchical_grn_moments()],
#'   [simulate_bulk_from_cells()],
#'   [single_cell_scenario_design()]
rescale_single_cell_moments <- function(
  moments,
  cell_size,
  reference = c("cellular_raw", "transcript_raw")
) {
  reference <- match.arg(reference)
  if (
    !is.list(moments) ||
      is.null(moments$mean_profiles) ||
      is.null(moments$covariance_matrices)
  ) {
    stop(
      "`moments` must contain `mean_profiles` and ",
      "`covariance_matrices`.",
      call. = FALSE
    )
  }
  mu <- as.matrix(moments$mean_profiles)
  sigma <- moments$covariance_matrices
  n_genes <- nrow(mu)
  n_types <- ncol(mu)
  if (
    length(dim(sigma)) != 3L ||
      !identical(dim(sigma)[1:2], c(n_genes, n_genes)) ||
      dim(sigma)[[3L]] != n_types
  ) {
    stop(
      "`covariance_matrices` must be a G x G x J array ",
      "matching `mean_profiles`.",
      call. = FALSE
    )
  }
  if (any(!is.finite(mu)) || any(!is.finite(sigma))) {
    stop("`moments` must be finite.", call. = FALSE)
  }
  if (any(mu < -1e-8)) {
    stop(
      "Means must be nonnegative. Call ",
      "simulate_hierarchical_grn_moments(..., nonnegative = TRUE).",
      call. = FALSE
    )
  }
  mu[mu < 0] <- 0
  column_size <- colSums(mu)
  if (any(column_size <= 0)) {
    stop("Each reference column must have a positive sum.", call. = FALSE)
  }

  cell_size <- as.numeric(cell_size)
  if (length(cell_size) == 1L) {
    cell_size <- rep(cell_size, n_types)
  }
  if (
    length(cell_size) != n_types ||
      any(!is.finite(cell_size)) ||
      any(cell_size <= 0)
  ) {
    stop(
      "`cell_size` must be a positive scalar or a length-J vector.",
      call. = FALSE
    )
  }
  if (!is.null(colnames(mu))) {
    names(cell_size) <- colnames(mu)
  }

  theta <- sweep(mu, 2L, column_size, "/")
  cellular_mu <- sweep(theta, 2L, cell_size, "*")
  dimnames(cellular_mu) <- dimnames(mu)
  cellular_sigma <- sigma
  dimnames(cellular_sigma) <- list(
    rownames(mu),
    rownames(mu),
    colnames(mu)
  )

  if (reference == "cellular_raw") {
    fit_mu <- cellular_mu
    fit_sigma <- cellular_sigma
    estimand <- "cellular"
  } else {
    fit_mu <- theta
    dimnames(fit_mu) <- dimnames(mu)
    fit_sigma <- cellular_sigma
    for (j in seq_len(n_types)) {
      fit_sigma[,, j] <- cellular_sigma[,, j] / cell_size[[j]]
    }
    estimand <- "transcript"
  }

  list(
    mean_signature_matrix = fit_mu,
    Sigma = fit_sigma,
    generative_mean = cellular_mu,
    generative_sigma = cellular_sigma,
    cell_size = cell_size,
    composition = theta,
    reference = reference,
    estimand = estimand
  )
}

#' Four single-cell simulation scenarios
#'
#' @description
#' The scenario grid used to check the single-cell fit. Two
#' compositions, equal and imbalanced, crossed with equal and unequal
#' cell sizes. Every scenario is drawn at \eqn{G=10} genes,
#' \eqn{n_i=10^{4}} cells and platform scales \eqn{a_i\in\{1, 2\}}
#' (same capture efficiency, or twice the bulk capture of the
#' single-cell platform).
#'
#' @return A list with `n_genes`, `n_cells`, `platform_scale` and
#'   `scenarios` (four named lists with `proportions` and `cell_size`).
#'
#' @examples
#' names(single_cell_scenario_design()$scenarios)
#' @export
#' @seealso [rescale_single_cell_moments()],
#'   [simulate_bulk_from_cells()]
single_cell_scenario_design <- function() {
  list(
    n_genes = 10L,
    n_cells = 10000,
    platform_scale = c(1, 2),
    scenarios = list(
      balanced_equal_size = list(
        proportions = c(0.5, 0.5),
        cell_size = c(1, 1)
      ),
      balanced_varying_size = list(
        proportions = c(0.5, 0.5),
        cell_size = c(1, 5)
      ),
      imbalanced_equal_size = list(
        proportions = c(0.8, 0.2),
        cell_size = c(1, 1)
      ),
      imbalanced_varying_size = list(
        proportions = c(0.8, 0.2),
        cell_size = c(1, 5)
      )
    )
  )
}
