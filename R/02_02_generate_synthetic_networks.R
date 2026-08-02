# ---- simulate_synthetic_GRNs.R --------------------------------------------
##
## Generative model for Gaussian-based gene regulatory networks with
## AutoGeneS-inspired mean profiles and graph-constrained precision
## matrices.

#' Score mean profiles with AutoGeneS-style objectives
#'
#' Computes the two objectives used by AutoGeneS for signature quality:
#' mean absolute pairwise cosine similarity between cell-type columns of
#' \eqn{\boldsymbol{\mu}} (to minimise) and the sum of pairwise Euclidean
#' distances between those columns (to maximise).
#'
#' @param mean_signature_matrix Numeric matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (columns = cell types).
#'
#' @return A named list with \code{mean_abs_cosine} and
#'   \code{sum_euclidean_distance}.
#'
#' @keywords internal
#' @noRd
compute_mean_profile_objectives <- function(mean_signature_matrix) {
  n_celltypes <- ncol(mean_signature_matrix)
  if (n_celltypes < 2L) {
    stop("`mean_signature_matrix` must have at least two columns.")
  }

  mean_abs_cosine <- 0
  sum_euclidean_distance <- 0
  n_pairs <- 0L

  for (j in seq_len(n_celltypes - 1L)) {
    for (k in (j + 1L):n_celltypes) {
      mu_j <- mean_signature_matrix[, j]
      mu_k <- mean_signature_matrix[, k]
      norms <- sqrt(sum(mu_j^2) * sum(mu_k^2))
      cosine_jk <- if (norms < .Machine$double.eps) {
        0
      } else {
        sum(mu_j * mu_k) / norms
      }
      mean_abs_cosine <- mean_abs_cosine + abs(cosine_jk)
      sum_euclidean_distance <- sum_euclidean_distance +
        sqrt(sum((mu_j - mu_k)^2))
      n_pairs <- n_pairs + 1L
    }
  }

  list(
    mean_abs_cosine = mean_abs_cosine / n_pairs,
    sum_euclidean_distance = sum_euclidean_distance
  )
}


#' Generate mean profiles with controlled cosine similarity
#'
#' Builds positive \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} by
#' blending cell-type-private marker blocks with a shared background
#' direction,
#' \deqn{
#'   \boldsymbol{\mu}_{\cdot j}
#'   = s\,
#'   \frac{
#'     \sqrt{\rho}\,\boldsymbol{u}
#'     +\sqrt{1-\rho}\,\boldsymbol{v}_{j}
#'   }{
#'     \bigl\|
#'       \sqrt{\rho}\,\boldsymbol{u}
#'       +\sqrt{1-\rho}\,\boldsymbol{v}_{j}
#'     \bigr\|_2
#'   },
#' }
#' where \eqn{\boldsymbol{v}_{j}} is high on the gene block of type \eqn{j}
#' and low elsewhere, \eqn{\boldsymbol{u}} is the shared (uniform) direction,
#' \eqn{\rho\in[0,1]} is \code{target_cosine}, and \eqn{s>0} is
#' \code{mean_scale}. Small \eqn{\rho} yields nearly complementary centroids;
#' large \eqn{\rho} recovers collinear profiles; small \eqn{s} recovers the
#' noise-sensitive low-distance regime.
#'
#' @param n_genes Integer \eqn{G}; must be at least \code{n_celltypes}.
#' @param n_celltypes Integer \eqn{J\ge 2}.
#' @param mean_scale Positive scalar \eqn{s} controlling centroid norms
#'   (hence Euclidean separation).
#' @param target_cosine Numeric in \eqn{[0,1]}, collinearity dial
#'   \eqn{\rho}.
#' @param gene_names Optional character vector of length \eqn{G}.
#' @param celltype_names Optional character vector of length \eqn{J}.
#' @param background_level Positive background level on non-marker genes
#'   within each private direction (default \code{0.05}).
#'
#' @return Numeric matrix \eqn{\boldsymbol{\mu}} with dimensions
#'   \eqn{G\times J}.
#'
#' @keywords internal
#' @noRd
generate_mean_signature_matrix <- function(
  n_genes,
  n_celltypes,
  mean_scale,
  target_cosine = 0,
  gene_names = NULL,
  celltype_names = NULL,
  background_level = 0.05
) {
  if (n_genes < n_celltypes) {
    stop("`n_genes` must be at least `n_celltypes`.")
  }
  if (n_celltypes < 2L) {
    stop("`n_celltypes` must be at least 2.")
  }
  if (mean_scale <= 0) {
    stop("`mean_scale` must be positive.")
  }
  if (target_cosine < 0 || target_cosine > 1) {
    stop("`target_cosine` must lie in [0, 1].")
  }
  if (background_level <= 0) {
    stop("`background_level` must be positive.")
  }

  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(n_genes))
  }
  if (is.null(celltype_names)) {
    celltype_names <- paste0("celltype_", seq_len(n_celltypes))
  }

  shared_direction <- rep(1 / sqrt(n_genes), n_genes)

  block_size <- n_genes %/% n_celltypes
  remainder <- n_genes - block_size * n_celltypes
  private_directions <- matrix(
    background_level,
    nrow = n_genes,
    ncol = n_celltypes
  )
  start_idx <- 1L
  for (j in seq_len(n_celltypes)) {
    len_j <- block_size + if (j <= remainder) 1L else 0L
    idx <- start_idx:(start_idx + len_j - 1L)
    private_directions[idx, j] <- 1
    private_directions[, j] <- private_directions[, j] /
      sqrt(sum(private_directions[, j]^2))
    start_idx <- start_idx + len_j
  }

  sqrt_rho <- sqrt(target_cosine)
  sqrt_one_minus_rho <- sqrt(1 - target_cosine)

  mean_signature_matrix <- matrix(
    NA_real_,
    nrow = n_genes,
    ncol = n_celltypes,
    dimnames = list(gene_names, celltype_names)
  )
  for (j in seq_len(n_celltypes)) {
    direction <- sqrt_rho *
      shared_direction +
      sqrt_one_minus_rho * private_directions[, j]
    direction <- direction / sqrt(sum(direction^2))
    mean_signature_matrix[, j] <- mean_scale * direction
  }

  mean_signature_matrix
}


#' Sample a symmetric adjacency matrix from a random graph model
#'
#' Wraps \pkg{igraph} generators for either a Barabási–Albert
#' (preferential-attachment / power-law) graph or a stochastic block
#' model.
#'
#' @param n_genes Integer, number of nodes (genes).
#' @param graph_model Character, one of \code{"power_law"} or
#'   \code{"stochastic_block_model"}.
#' @param graph_params Named list of model-specific parameters (see
#'   \code{\link{simulate_hierarchical_grn_moments}} for details).
#'
#' @return A symmetric integer matrix of dimension \eqn{G\times G} with
#'   zero diagonal.
#'
#' @keywords internal
#' @noRd
generate_random_network_skeleton <- function(
  n_genes,
  graph_model,
  graph_params
) {
  graph <- switch(
    graph_model,

    power_law = {
      power <- if (is.null(graph_params$power)) {
        1
      } else {
        graph_params$power
      }
      edges_per_node <- if (is.null(graph_params$edges_per_node)) {
        1L
      } else {
        as.integer(graph_params$edges_per_node)
      }
      igraph::sample_pa(
        n_genes,
        power = power,
        m = edges_per_node,
        directed = FALSE
      )
    },

    stochastic_block_model = {
      block_prob <- if (is.null(graph_params$block_prob)) {
        c(0.5, 0.25, 0.25)
      } else {
        graph_params$block_prob
      }
      p_within <- if (is.null(graph_params$p_within)) {
        0.25
      } else {
        graph_params$p_within
      }
      p_between <- if (is.null(graph_params$p_between)) {
        0.01
      } else {
        graph_params$p_between
      }

      n_blocks <- length(block_prob)
      pref_matrix <- matrix(p_between, n_blocks, n_blocks)
      diag(pref_matrix) <- p_within

      block_sizes <- as.integer(
        c(stats::rmultinom(1L, n_genes, block_prob))
      )
      block_sizes <- pmax(block_sizes, 1L)
      size_diff <- n_genes - sum(block_sizes)
      if (size_diff != 0L) {
        idx <- which.max(block_sizes)
        block_sizes[idx] <- block_sizes[idx] + size_diff
      }

      igraph::sample_sbm(
        n_genes,
        pref.matrix = pref_matrix,
        block.sizes = block_sizes
      )
    }
  )

  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(igraph::as_undirected(graph))
  )
  diag(adjacency) <- 0L
  adjacency
}


#' Build a normalised precision matrix from a binary adjacency matrix
#'
#' Applies an affine transformation to the adjacency matrix so that the
#' resulting precision matrix is symmetric positive-definite:
#' \deqn{
#'   \tilde{\boldsymbol{\Omega}} = \boldsymbol{A} \odot v, \quad
#'   \boldsymbol{\Omega}
#'   = \tilde{\boldsymbol{\Omega}}
#'   + \bigl(\lvert\lambda_{\min}(\tilde{\boldsymbol{\Omega}})\rvert + u\bigr)
#'     \mathbf{I}
#' }
#'
#' @param adjacency_matrix Symmetric binary matrix.
#' @param precision_shift Positive numeric (\eqn{u}), additive diagonal
#'   shift above the absolute smallest eigenvalue.
#' @param precision_scale Positive numeric (\eqn{v}), multiplicative
#'   scale for off-diagonal entries.
#'
#' @return A symmetric positive-definite matrix of the same dimension.
#'
#' @keywords internal
#' @noRd
build_normalised_precision <- function(
  adjacency_matrix,
  precision_shift,
  precision_scale
) {
  omega <- adjacency_matrix * precision_scale
  eigenvalues <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
  diag(omega) <- abs(min(eigenvalues)) + precision_shift
  omega
}


#' Replicate a shared precision as cell-type covariances
#'
#' Inverts \eqn{\boldsymbol{\Omega}} to
#' \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Omega}^{-1}} and stacks the same
#' covariance for each of the \eqn{J} cell types, matching the article's
#' shared graph-constrained second-order structure.
#'
#' @param normalised_precision Symmetric positive-definite
#'   \eqn{\boldsymbol{\Omega}\in\mathcal{M}_{G\times G}}.
#' @param celltype_names Character vector of length \eqn{J}.
#' @param gene_names Character vector of length \eqn{G}.
#'
#' @return Numeric array of dimension \eqn{G\times G\times J}.
#'
#' @keywords internal
#' @noRd
build_shared_covariance_array <- function(
  normalised_precision,
  celltype_names,
  gene_names
) {
  n_genes <- nrow(normalised_precision)
  n_celltypes <- length(celltype_names)
  sigma <- solve(normalised_precision)
  sigma <- (sigma + t(sigma)) / 2

  covariance_array <- array(
    NA_real_,
    dim = c(n_genes, n_genes, n_celltypes),
    dimnames = list(gene_names, gene_names, celltype_names)
  )
  for (j in seq_len(n_celltypes)) {
    covariance_array[,, j] <- sigma
  }
  covariance_array
}


#' Simulate GRN first- and second-order moments
#'
#' @description
#' Builds a mean matrix \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}
#' and a shared covariance array \eqn{(\boldsymbol{\Sigma}_j)_{j}} under a
#' graph-constrained precision model. Means follow the AutoGeneS-inspired
#' construction of `generate_mean_signature_matrix()` with target pairwise
#' cosine similarity \eqn{\rho} and column scale \eqn{s}, so that one can
#' dial collinearity (minimise cosine) against centroid separation
#' (maximise Euclidean distance). The adjacency
#' \eqn{\boldsymbol{A}} is drawn from a random-graph model; the precision
#' \eqn{\boldsymbol{\Omega}} is obtained via the affine spectral shift
#' of `build_normalised_precision()`, and each cell type shares
#' \eqn{\boldsymbol{\Sigma}_j=\boldsymbol{\Omega}^{-1}}.
#'
#' @param n_genes Integer; number of genes \eqn{G}.
#' @param n_celltypes Integer; number of cell types \eqn{J} (default 2).
#' @param mean_scale Positive scalar \eqn{s} for centroid norms.
#' @param target_cosine Numeric in \eqn{[0,1]}; target pairwise cosine
#'   similarity between columns of \eqn{\boldsymbol{\mu}}.
#' @param precision_shift,precision_scale Diagonal shift \eqn{u} and
#'   off-diagonal scale \eqn{v} used to build \eqn{\boldsymbol{\Omega}}.
#' @param graph_model `"power_law"` or `"stochastic_block_model"`.
#' @param graph_params Named list of generator parameters:
#'   `power` / `edges_per_node` (power-law) or
#'   `block_prob` / `p_within` / `p_between` (SBM).
#'
#' @return Named list with:
#' * `mean_profiles`: matrix \eqn{\boldsymbol{\mu}};
#' * `covariance_matrices`: array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}};
#' * `graph_structure`: `adjacency_matrix` and `normalised_precision`
#'   \eqn{\boldsymbol{\Omega}};
#' * `objectives`: `mean_abs_cosine` and `sum_euclidean_distance`.
#'
#' @examples
#' set.seed(42)
#' moments <- simulate_hierarchical_grn_moments(
#'   n_genes = 40L,
#'   n_celltypes = 3L,
#'   mean_scale = 10,
#'   target_cosine = 0.1,
#'   precision_shift = 0.1,
#'   precision_scale = 0.3,
#'   graph_model = "power_law",
#'   graph_params = list(power = 1, edges_per_node = 2)
#' )
#' str(moments, max.level = 2)
#'
#' ## Verify positive-definiteness of the shared covariance
#' eigen_vals <- eigen(
#'   moments$covariance_matrices[, , 1],
#'   only.values = TRUE
#' )$values
#' stopifnot(all(eigen_vals > 0))
#'
#' @importFrom igraph sample_pa sample_sbm as_adjacency_matrix
#'   as_undirected
#' @importFrom stats rmultinom
#'
#' @export
simulate_hierarchical_grn_moments <- function(
  n_genes,
  n_celltypes = 2L,
  mean_scale = 10,
  target_cosine = 0,
  precision_shift,
  precision_scale,
  graph_model = c("power_law", "stochastic_block_model"),
  graph_params = list()
) {
  graph_model <- match.arg(graph_model)
  n_genes <- as.integer(n_genes)
  n_celltypes <- as.integer(n_celltypes)

  if (n_genes < n_celltypes) {
    stop("`n_genes` must be at least `n_celltypes`.")
  }
  if (precision_shift <= 0 || precision_scale <= 0) {
    stop("`precision_shift` and `precision_scale` must be positive.")
  }

  gene_names <- paste0("gene_", seq_len(n_genes))
  celltype_names <- paste0("celltype_", seq_len(n_celltypes))

  ## --- step 1: AutoGeneS-inspired mean profiles ------------------------
  mean_profiles <- generate_mean_signature_matrix(
    n_genes = n_genes,
    n_celltypes = n_celltypes,
    mean_scale = mean_scale,
    target_cosine = target_cosine,
    gene_names = gene_names,
    celltype_names = celltype_names
  )
  objectives <- compute_mean_profile_objectives(mean_profiles)

  ## --- step 2: random-graph adjacency ----------------------------------
  adjacency_matrix <- generate_random_network_skeleton(
    n_genes = n_genes,
    graph_model = graph_model,
    graph_params = graph_params
  )
  rownames(adjacency_matrix) <- gene_names
  colnames(adjacency_matrix) <- gene_names

  ## --- step 3: normalised precision and shared covariances -------------
  normalised_precision <- build_normalised_precision(
    adjacency_matrix = adjacency_matrix,
    precision_shift = precision_shift,
    precision_scale = precision_scale
  )
  dimnames(normalised_precision) <- list(gene_names, gene_names)

  covariance_matrices <- build_shared_covariance_array(
    normalised_precision = normalised_precision,
    celltype_names = celltype_names,
    gene_names = gene_names
  )

  list(
    mean_profiles = mean_profiles,
    covariance_matrices = covariance_matrices,
    graph_structure = list(
      adjacency_matrix = adjacency_matrix,
      normalised_precision = normalised_precision
    ),
    objectives = objectives
  )
}
