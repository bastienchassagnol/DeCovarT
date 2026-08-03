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


#' Generate mean profiles with a target pairwise cosine
#'
#' Builds positive
#' \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} by blending a shared
#' unit direction \eqn{\boldsymbol{u}} with cell-type-private nearly
#' orthogonal marker directions \eqn{\boldsymbol{v}_{j}}:
#' \deqn{
#'   \tilde{\boldsymbol{\mu}}_{\cdot j}
#'   =
#'   \sqrt{\rho}\,\boldsymbol{u}
#'   +\sqrt{1-\rho}\,\boldsymbol{v}_{j},
#'   \qquad
#'   \boldsymbol{\mu}_{\cdot j}
#'   =
#'   s\,
#'   \frac{\tilde{\boldsymbol{\mu}}_{\cdot j}}{
#'     \|\tilde{\boldsymbol{\mu}}_{\cdot j}\|_2
#'   }.
#' }
#' Here \eqn{\rho\in[0,1]} is \code{target_cosine} and \eqn{s>0} is
#' \code{mean_scale}. The private vectors \eqn{\boldsymbol{v}_{j}} are
#' supported on a partition of the \eqn{G} genes (high on the block of
#' type \eqn{j}, small background elsewhere) and then
#' \eqn{\ell_2}-normalised, so
#' \eqn{\boldsymbol{v}_{j}^{\mathsf{T}}\boldsymbol{v}_{k}\approx 0} for
#' \eqn{j\neq k}. With a shared unit \eqn{\boldsymbol{u}},
#' \deqn{
#'   \tilde{\boldsymbol{\mu}}_{\cdot j}^{\mathsf{T}}
#'   \tilde{\boldsymbol{\mu}}_{\cdot k}
#'   \approx
#'   \rho
#'   \qquad (j\neq k),
#' }
#' and after column normalisation the pairwise cosines of
#' \eqn{\boldsymbol{\mu}} therefore track \eqn{\rho} (exactly in the
#' ideal orthonormal limit; approximately with a small background
#' level). The global scale \eqn{s} sets column norms (and hence
#' Euclidean separation) without changing angles: for fixed
#' \eqn{\rho},
#' \eqn{\|\boldsymbol{\mu}_{\cdot j}-\boldsymbol{\mu}_{\cdot k}\|_2
#' \propto s}. Prefer dialling \eqn{\rho} when second-order precision
#' weights already control interaction strength; keep \eqn{s} fixed
#' across scenarios that compare mean collinearity alone
#' \insertCite{alieeAutoGeneSAutomaticGene2021}{DeCovarT}.
#'
#' @param n_genes Integer \eqn{G}; must be at least \code{n_celltypes}.
#' @param n_celltypes Integer \eqn{J\ge 2}.
#' @param mean_scale Positive scalar \eqn{s} (centroid norms). Hold
#'   fixed when studying cosine / collinearity alone.
#' @param target_cosine Numeric in \eqn{[0,1]}, the collinearity dial
#'   \eqn{\rho}.
#' @param gene_names Optional character vector of length \eqn{G}.
#' @param celltype_names Optional character vector of length \eqn{J}.
#' @param background_level Positive background on non-marker genes in
#'   each private direction (default \code{0.05}).
#'
#' @return Numeric matrix \eqn{\boldsymbol{\mu}} with dimensions
#'   \eqn{G\times J}.
#'
#' @keywords internal
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


#' Watts–Strogatz small-world adjacency (undirected)
#'
#' Ring lattice with each node linked to \code{nei} neighbours on each
#' side, then undirected edge rewiring with probability \code{p}.
#' Used because \pkg{huge} does not expose a small-world generator.
#'
#' @keywords internal
.sample_small_world_adjacency <- function(n_genes, nei = 2L, p = 0.05) {
  n_genes <- as.integer(n_genes)
  nei <- as.integer(nei)
  if (nei < 1L || 2L * nei >= n_genes) {
    stop("`nei` must satisfy 1 <= nei < n_genes / 2.")
  }
  if (p < 0 || p > 1) {
    stop("`p` must lie in [0, 1].")
  }

  adjacency <- matrix(0L, n_genes, n_genes)
  for (i in seq_len(n_genes)) {
    for (d in seq_len(nei)) {
      j <- ((i - 1L + d) %% n_genes) + 1L
      adjacency[i, j] <- 1L
      adjacency[j, i] <- 1L
    }
  }

  upper <- which(upper.tri(adjacency) & adjacency == 1L, arr.ind = TRUE)
  if (nrow(upper) == 0L) {
    return(adjacency)
  }

  for (e in seq_len(nrow(upper))) {
    if (stats::runif(1L) > p) {
      next
    }
    i <- upper[e, 1L]
    j <- upper[e, 2L]
    candidates <- which(adjacency[i, ] == 0L & seq_len(n_genes) != i)
    if (length(candidates) == 0L) {
      next
    }
    k <- candidates[[sample.int(length(candidates), 1L)]]
    adjacency[i, j] <- 0L
    adjacency[j, i] <- 0L
    adjacency[i, k] <- 1L
    adjacency[k, i] <- 1L
  }

  diag(adjacency) <- 0L
  adjacency
}


#' Sample a symmetric adjacency from a random-graph family
#'
#' Preferential-attachment (scale-free) and cluster / stochastic-block
#' skeletons are drawn with \code{huge::huge.generator()}
#' \insertCite{zhangSILGGMExtensivePackage2018}{DeCovarT}. Small-world
#' graphs use an internal Watts–Strogatz construction because \pkg{huge}
#' does not expose that family.
#'
#' @param n_genes Integer \eqn{G}, number of nodes (genes).
#' @param graph_model One of \code{"scale_free"},
#'   \code{"stochastic_block_model"}, \code{"small_world"}.
#' @param graph_params Named list of generator parameters:
#'   \describe{
#'     \item{\code{scale_free}}{\code{verbose} (logical)}
#'     \item{\code{stochastic_block_model}}{\code{n_blocks} (passed as
#'       \code{g} to \code{huge.generator}), \code{verbose}}
#'     \item{\code{small_world}}{\code{nei}, \code{p} (rewiring
#'       probability)}
#'   }
#'
#' @return Symmetric integer matrix \eqn{G\times G} with zero diagonal.
#'
#' @keywords internal
generate_random_network_skeleton <- function(
  n_genes,
  graph_model,
  graph_params = list()
) {
  n_genes <- as.integer(n_genes)
  graph_model <- match.arg(
    graph_model,
    c("scale_free", "stochastic_block_model", "small_world")
  )

  adjacency <- switch(
    graph_model,

    scale_free = {
      verbose <- isTRUE(graph_params$verbose)
      sim <- huge::huge.generator(
        n = max(2L, n_genes),
        d = n_genes,
        graph = "scale-free",
        vis = FALSE,
        verbose = verbose
      )
      as.matrix(sim$theta)
    },

    stochastic_block_model = {
      n_blocks <- if (is.null(graph_params$n_blocks)) {
        3L
      } else {
        as.integer(graph_params$n_blocks)
      }
      verbose <- isTRUE(graph_params$verbose)
      sim <- huge::huge.generator(
        n = max(2L, n_genes),
        d = n_genes,
        graph = "cluster",
        g = n_blocks,
        vis = FALSE,
        verbose = verbose
      )
      as.matrix(sim$theta)
    },

    small_world = {
      nei <- if (is.null(graph_params$nei)) {
        2L
      } else {
        as.integer(graph_params$nei)
      }
      p_rewire <- if (is.null(graph_params$p)) {
        0.05
      } else {
        graph_params$p
      }
      .sample_small_world_adjacency(
        n_genes = n_genes,
        nei = nei,
        p = p_rewire
      )
    }
  )

  storage.mode(adjacency) <- "integer"
  adjacency <- (adjacency + t(adjacency))
  adjacency[adjacency > 0L] <- 1L
  diag(adjacency) <- 0L
  adjacency
}


#' Assign i.i.d. signed weights on a binary undirected adjacency
#'
#' For each undirected edge, draws an inhibitory versus activatory sign
#' and a common magnitude. By convention a **positive** precision
#' off-diagonal is an **inhibitory** partial correlation
#' (\eqn{\rho_{jk\mid -\{j,k\}}=-\Omega_{jk}/
#' \sqrt{\Omega_{jj}\Omega_{kk}}}). Thus
#' \code{prop_inhibitory} is the target fraction of edges with
#' \eqn{\Omega_{jk}>0} (negative partial correlation). The complementary
#' fraction is activatory (\eqn{\Omega_{jk}<0}).
#'
#' @param adjacency_matrix Symmetric binary matrix.
#' @param prop_inhibitory Numeric in \eqn{[0,1]}; expected fraction of
#'   inhibitory precision edges.
#' @param weight_magnitude Positive magnitude \eqn{v} for every signed
#'   off-diagonal.
#'
#' @return Symmetric numeric matrix with the same support as
#'   \code{adjacency_matrix} and zero diagonal.
#'
#' @keywords internal
assign_iid_signed_weights <- function(
  adjacency_matrix,
  prop_inhibitory = 0.5,
  weight_magnitude = 0.3
) {
  A <- as.matrix(adjacency_matrix)
  if (nrow(A) != ncol(A)) {
    stop("`adjacency_matrix` must be square.")
  }
  if (prop_inhibitory < 0 || prop_inhibitory > 1) {
    stop("`prop_inhibitory` must lie in [0, 1].")
  }
  if (weight_magnitude <= 0) {
    stop("`weight_magnitude` must be positive.")
  }

  W <- matrix(0, nrow(A), ncol(A))
  upper <- which(upper.tri(A) & A != 0, arr.ind = TRUE)
  if (nrow(upper) == 0L) {
    dimnames(W) <- dimnames(A)
    return(W)
  }

  n_edge <- nrow(upper)
  n_inhib <- as.integer(round(prop_inhibitory * n_edge))
  signs <- rep(-1, n_edge) # activatory precision weight by default
  if (n_inhib > 0L) {
    inhib_idx <- sample.int(n_edge, n_inhib)
    signs[inhib_idx] <- 1
  }
  values <- signs * weight_magnitude
  W[upper] <- values
  W[cbind(upper[, 2L], upper[, 1L])] <- values
  dimnames(W) <- dimnames(A)
  W
}


#' Complete a weighted adjacency to an SPD precision
#'
#' Applies a uniform spectral shift that preserves off-diagonal signs and
#' support:
#' \deqn{
#'   \boldsymbol{\Omega}
#'   =
#'   \boldsymbol{W}
#'   +
#'   \bigl(\lvert\lambda_{\min}(\boldsymbol{W})\rvert + u\bigr)
#'   \mathbf{I}.
#' }
#'
#' @param weighted_adjacency Symmetric numeric matrix (zero diagonal).
#' @param precision_shift Positive diagonal cushion \eqn{u}.
#'
#' @return Symmetric positive-definite matrix.
#'
#' @keywords internal
build_normalised_precision <- function(
  weighted_adjacency,
  precision_shift
) {
  if (precision_shift <= 0) {
    stop("`precision_shift` must be positive.")
  }
  omega <- as.matrix(weighted_adjacency)
  diag(omega) <- 0
  eigenvalues <- eigen(omega, symmetric = TRUE, only.values = TRUE)$values
  diag(omega) <- abs(min(eigenvalues)) + precision_shift
  omega
}


#' Replicate a shared precision as cell-type covariances
#'
#' Inverts \eqn{\boldsymbol{\Omega}} to
#' \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Omega}^{-1}} and stacks the same
#' covariance for each of the \eqn{J} cell types.
#'
#' @param normalised_precision Symmetric positive-definite
#'   \eqn{\boldsymbol{\Omega}\in\mathcal{M}_{G\times G}}.
#' @param celltype_names Character vector of length \eqn{J}.
#' @param gene_names Character vector of length \eqn{G}.
#'
#' @return Numeric array of dimension \eqn{G\times G\times J}.
#'
#' @keywords internal
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
#' cosine construction of `generate_mean_signature_matrix()` with target
#' pairwise cosine \eqn{\rho}. The adjacency \eqn{\boldsymbol{A}} is drawn
#' from a random-graph model (\pkg{huge} for scale-free and cluster /
#' SBM; internal Watts–Strogatz for small-world); i.i.d. signed weights
#' with inhibitory fraction \code{prop_inhibitory} form
#' \eqn{\boldsymbol{W}}; the precision is completed by a spectral shift;
#' each cell type shares
#' \eqn{\boldsymbol{\Sigma}_j=\boldsymbol{\Omega}^{-1}}.
#'
#' @param n_genes Integer; number of genes \eqn{G}.
#' @param n_celltypes Integer; number of cell types \eqn{J} (default 2).
#' @param mean_scale Positive scalar \eqn{s} for centroid norms.
#' @param target_cosine Numeric in \eqn{[0,1]}; target pairwise cosine
#'   similarity between columns of \eqn{\boldsymbol{\mu}}.
#' @param precision_shift Diagonal cushion \eqn{u} for the spectral
#'   shift.
#' @param precision_scale Positive magnitude \eqn{v} of signed
#'   off-diagonal precision weights.
#' @param prop_inhibitory Numeric in \eqn{[0,1]}; fraction of edges with
#'   positive precision weight (inhibitory partial correlation). Default
#'   \code{0.5} balances inhibitory and activatory edges.
#' @param graph_model One of \code{"scale_free"},
#'   \code{"stochastic_block_model"}, \code{"small_world"}.
#' @param graph_params Named list of generator parameters (see
#'   \code{generate_random_network_skeleton()}).
#'
#' @return Named list with:
#' * `mean_profiles`: matrix \eqn{\boldsymbol{\mu}};
#' * `covariance_matrices`: array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}};
#' * `graph_structure`: `adjacency_matrix`, `weighted_adjacency`, and
#'   `normalised_precision` \eqn{\boldsymbol{\Omega}};
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
#'   prop_inhibitory = 0.5,
#'   graph_model = "scale_free"
#' )
#' str(moments, max.level = 2)
#'
#' @importFrom stats runif
#'
#' @export
simulate_hierarchical_grn_moments <- function(
  n_genes,
  n_celltypes = 2L,
  mean_scale = 10,
  target_cosine = 0,
  precision_shift,
  precision_scale,
  prop_inhibitory = 0.5,
  graph_model = c(
    "scale_free",
    "stochastic_block_model",
    "small_world"
  ),
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

  ## --- step 3: i.i.d. signs then spectral SPD completion ---------------
  weighted_adjacency <- assign_iid_signed_weights(
    adjacency_matrix = adjacency_matrix,
    prop_inhibitory = prop_inhibitory,
    weight_magnitude = precision_scale
  )
  normalised_precision <- build_normalised_precision(
    weighted_adjacency = weighted_adjacency,
    precision_shift = precision_shift
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
      weighted_adjacency = weighted_adjacency,
      normalised_precision = normalised_precision
    ),
    objectives = objectives
  )
}
