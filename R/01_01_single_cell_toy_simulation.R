# Single-cell bulk: a sum of independent Gaussian cells.
# The population generator simulate_bulk_mixture() (p^2 weights) is in
# R/01_01_population_toy_simulation.R and is not used here. The fit needs
# G > J genes to identify the two technical scales.

#' Simulate a bulk sample as a sum of single cells
#'
#' @description
#' Draws cell counts
#' \eqn{m_{\cdot i}\sim Multinomial(n_i,p_{\cdot i})},
#' then the exact Gaussian sum of those cells,
#' \eqn{T_{ji}\sim N(m_{ji}\mu_{\cdot j},m_{ji}\Sigma_{W,j})},
#' and scales the total by the platform factor:
#' \eqn{\boldsymbol{y}_{\cdot i}=a_i\sum_j\boldsymbol{T}_{ji}}.
#' A sum of \eqn{m} independent Gaussian cells is Gaussian, so there is
#' no per-cell loop.
#'
#' This is the single-cell mixture (covariance linear in the cell
#' counts). [simulate_bulk_mixture()] remains the population-profile
#' convolution, whose covariance is weighted by \eqn{p_j^2}.
#'
#' Pass the cellular moments (`generative_mean`, `generative_sigma`
#' from [rescale_single_cell_moments()]), not a library-normalised
#' reference. The returned estimand is the realised cell fraction
#' \eqn{\boldsymbol{m}_{\cdot i}/n_i}.
#'
#' @param mean_signature_matrix Per-cell means
#'   \eqn{\tilde{\boldsymbol{\mu}}} (\eqn{G\times J}).
#' @param Sigma Per-cell covariances
#'   \eqn{\boldsymbol{\Sigma}_{W,j}} (\eqn{G\times G\times J}).
#' @param p Cell-fraction vector of length \eqn{J}, or a \eqn{J\times N}
#'   matrix. Each column must sum to 1.
#' @param n_cells Cells per bulk sample, \eqn{n_i}. Recycled to \eqn{N}.
#'   The multinomial size is the rounded value; the fitted
#'   \eqn{\hat n} need not be an integer.
#' @param platform_scale Capture scale \eqn{a_i}. Use `1` when the bulk
#'   and the single-cell platform share their scale, and `2` when the
#'   bulk captures twice as much transcript per cell. Recycled to
#'   \eqn{N}.
#' @param n_samples Number of bulk samples when `p` is a vector.
#'   Ignored when `p` already has one column per sample, unless it
#'   matches that number.
#' @param seed Optional seed, passed to [set.seed()] before the draws.
#'
#' @return A list with `Y` (\eqn{G\times N}), nominal `p`, realised
#'   `p_realised` (\eqn{\boldsymbol{m}/n}), `counts` (\eqn{\boldsymbol{m}}),
#'   transcript fractions `q`, `cell_size` \eqn{S_j=\mathbf{1}^{\top}
#'   \tilde{\boldsymbol{\mu}}_{\cdot j}}, `platform_scale` \eqn{a},
#'   `n_cells` \eqn{n}, `mean_scale` \eqn{c=an} and `variance_scale`
#'   \eqn{\phi=a^{2}n}.
#'
#' @examples
#' genes <- paste0("g", 1:4)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(
#'   c(2, 2, 0.2, 0.2, 0.2, 0.2, 3, 3),
#'   nrow = 4,
#'   dimnames = list(genes, cts)
#' )
#' sigma <- array(0, c(4, 4, 2), dimnames = list(genes, genes, cts))
#' sigma[,, 1] <- diag(0.05, 4)
#' sigma[,, 2] <- diag(0.05, 4)
#' set.seed(1)
#' sim <- simulate_bulk_from_cells(
#'   mu,
#'   sigma,
#'   p = c(0.8, 0.2),
#'   n_cells = 2000,
#'   platform_scale = 2,
#'   n_samples = 2
#' )
#' colSums(sim$counts)
#' @export
#' @seealso [simulate_bulk_mixture()],
#'   [rescale_single_cell_moments()],
#'   [deconvolute_single_cell()]
simulate_bulk_from_cells <- function(
  mean_signature_matrix,
  Sigma,
  p,
  n_cells,
  platform_scale = 1,
  n_samples = NULL,
  seed = NULL
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  mean_signature_matrix <- as.matrix(mean_signature_matrix)
  n_genes <- nrow(mean_signature_matrix)
  n_types <- ncol(mean_signature_matrix)
  if (n_genes <= n_types) {
    stop(
      "The two technical scales need more genes than cell types ",
      "(G > J); got G = ",
      n_genes,
      " and J = ",
      n_types,
      ".",
      call. = FALSE
    )
  }
  if (
    length(dim(Sigma)) != 3L ||
      !identical(dim(Sigma)[1:2], c(n_genes, n_genes)) ||
      dim(Sigma)[[3L]] != n_types
  ) {
    stop("`Sigma` must be a G x G x J array.", call. = FALSE)
  }
  if (is.null(rownames(mean_signature_matrix))) {
    rownames(mean_signature_matrix) <- paste0("gene_", seq_len(n_genes))
  }
  if (is.null(colnames(mean_signature_matrix))) {
    colnames(mean_signature_matrix) <- paste0(
      "celltype_",
      seq_len(n_types)
    )
  }
  gene_names <- rownames(mean_signature_matrix)
  cell_names <- colnames(mean_signature_matrix)
  sigma_genes <- dimnames(Sigma)[[1L]]
  sigma_cells <- dimnames(Sigma)[[3L]]
  if (!is.null(sigma_genes) && !identical(sigma_genes, gene_names)) {
    stop("Gene names differ between the mean and `Sigma`.", call. = FALSE)
  }
  if (!is.null(sigma_cells) && !identical(sigma_cells, cell_names)) {
    stop(
      "Cell-type names differ between the mean and `Sigma`.",
      call. = FALSE
    )
  }

  if (is.null(dim(p))) {
    p <- as.numeric(p)
    if (is.null(n_samples)) {
      n_samples <- 1L
    }
    n_samples <- as.integer(n_samples)
    p <- matrix(p, nrow = length(p), ncol = n_samples)
  } else {
    p <- as.matrix(p)
    if (is.null(n_samples)) {
      n_samples <- ncol(p)
    }
    n_samples <- as.integer(n_samples)
    if (n_samples != ncol(p)) {
      stop("`n_samples` does not match `ncol(p)`.", call. = FALSE)
    }
  }
  if (nrow(p) != n_types) {
    stop("`p` must have one row per cell type.", call. = FALSE)
  }
  if (any(!is.finite(p)) || any(p < -1e-8)) {
    stop("`p` must be a finite nonnegative composition.", call. = FALSE)
  }
  p[p < 0] <- 0
  column_total <- colSums(p)
  if (any(column_total <= 0) || any(abs(column_total - 1) > 1e-6)) {
    stop("Each column of `p` must sum to 1.", call. = FALSE)
  }
  p <- sweep(p, 2L, column_total, "/")
  rownames(p) <- cell_names

  n_draw <- .recycle_positive_samples(n_cells, n_samples, "n_cells")
  n_draw <- as.integer(round(n_draw))
  if (any(n_draw < 1L)) {
    stop("`n_cells` must be at least 1.", call. = FALSE)
  }
  platform_scale <- .recycle_positive_samples(
    platform_scale,
    n_samples,
    "platform_scale"
  )

  sample_names <- paste0("sample_", seq_len(n_samples))
  counts <- matrix(
    0L,
    nrow = n_types,
    ncol = n_samples,
    dimnames = list(cell_names, sample_names)
  )
  y <- matrix(
    0,
    nrow = n_genes,
    ncol = n_samples,
    dimnames = list(gene_names, sample_names)
  )
  # One Cholesky factor per type: count_j cells contribute
  # N(count_j mu_j, count_j Sigma_j), drawn as count_j mu_j + sqrt(count_j) z R.
  roots <- lapply(seq_len(n_types), function(j) .covariance_root(Sigma[,, j]))
  for (i in seq_len(n_samples)) {
    counts[, i] <- as.integer(stats::rmultinom(
      1L,
      n_draw[[i]],
      p[, i]
    ))
    total <- rep(0, n_genes)
    for (j in seq_len(n_types)) {
      count_j <- counts[[j, i]]
      if (count_j <= 0L) {
        next
      }
      draw <- count_j * mean_signature_matrix[, j] +
        sqrt(count_j) * drop(crossprod(roots[[j]], stats::rnorm(n_genes)))
      total <- total + as.numeric(draw)
    }
    y[, i] <- platform_scale[[i]] * total
  }

  p_realised <- sweep(counts, 2L, n_draw, "/")
  cell_size <- colSums(mean_signature_matrix)
  names(cell_size) <- cell_names
  size_mass <- cell_size * p_realised
  q <- sweep(size_mass, 2L, colSums(size_mass), "/")
  dimnames(q) <- dimnames(counts)
  dimnames(p_realised) <- dimnames(counts)

  list(
    Y = y,
    p = p,
    p_realised = p_realised,
    counts = counts,
    q = q,
    cell_size = cell_size,
    platform_scale = platform_scale,
    n_cells = n_draw,
    mean_scale = platform_scale * n_draw,
    variance_scale = platform_scale^2 * n_draw
  )
}

#' @noRd
.recycle_positive_samples <- function(x, n_samples, name) {
  x <- as.numeric(x)
  if (length(x) == 1L) {
    x <- rep(x, n_samples)
  }
  if (length(x) != n_samples || any(!is.finite(x)) || any(x <= 0)) {
    stop(
      "`",
      name,
      "` must be a positive scalar or a length-N vector.",
      call. = FALSE
    )
  }
  x
}
