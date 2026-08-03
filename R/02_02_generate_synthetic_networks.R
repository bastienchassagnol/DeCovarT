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
#' Builds
#' \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} by blending a shared
#' unit direction \eqn{\boldsymbol{u}} with cell-type-private orthogonal
#' marker directions \eqn{\boldsymbol{v}_{j}}:
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
#' \code{mean_scale} (default \code{10}, matching the nine-scenario
#' simulation grid). The private vectors \eqn{\boldsymbol{v}_{j}} are
#' indicator directions on a partition of the \eqn{G} genes (type
#' \eqn{j} only) and then \eqn{\ell_2}-normalised, so
#' \eqn{\boldsymbol{v}_{j}^{\mathsf{T}}\boldsymbol{v}_{k}=0} for
#' \eqn{j\neq k}. With a shared unit \eqn{\boldsymbol{u}},
#' \deqn{
#'   \tilde{\boldsymbol{\mu}}_{\cdot j}^{\mathsf{T}}
#'   \tilde{\boldsymbol{\mu}}_{\cdot k}
#'   =
#'   \rho
#'   +
#'   \sqrt{\rho(1-\rho)}\,
#'   \bigl(
#'     \boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}_{j}
#'     +
#'     \boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}_{k}
#'   \bigr)
#'   \qquad (j\neq k).
#' }
#' After column normalisation the pairwise cosines of
#' \eqn{\boldsymbol{\mu}} therefore track \eqn{\rho} closely when the
#' cross terms
#' \eqn{\boldsymbol{u}^{\mathsf{T}}\boldsymbol{v}_{j}} are small relative
#' to the leading \eqn{\rho} (many genes per block). The global scale
#' \eqn{s} sets column norms (and hence Euclidean separation) without
#' changing angles: for fixed \eqn{\rho},
#' \eqn{\|\boldsymbol{\mu}_{\cdot j}-\boldsymbol{\mu}_{\cdot k}\|_2
#' \propto s}. Prefer dialling \eqn{\rho} when second-order precision
#' weights already control interaction strength; keep \eqn{s} fixed
#' across scenarios that compare mean collinearity alone
#' \insertCite{alieeAutoGeneSAutomaticGene2021}{DeCovarT}.
#'
#' @param n_genes Integer \eqn{G}; must be at least \code{n_celltypes}.
#' @param n_celltypes Integer \eqn{J\ge 2}.
#' @param mean_scale Positive scalar \eqn{s} (centroid norms). Default
#'   \code{10}, as in the nine factorial scenarios. Hold fixed when
#'   studying cosine / collinearity alone.
#' @param target_cosine Numeric in \eqn{[0,1]}, the collinearity dial
#'   \eqn{\rho}.
#' @param gene_names Optional character vector of length \eqn{G}.
#' @param celltype_names Optional character vector of length \eqn{J}.
#'
#' @return Numeric matrix \eqn{\boldsymbol{\mu}} with dimensions
#'   \eqn{G\times J}.
#'
#' @keywords internal
generate_mean_signature_matrix <- function(
  n_genes,
  n_celltypes,
  mean_scale = 10,
  target_cosine = 0,
  gene_names = NULL,
  celltype_names = NULL
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
    0,
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


#' Sample a symmetric adjacency from a random-graph family
#'
#' Draws undirected skeletons with \pkg{igraph}: Barabási–Albert
#' preferential attachment (\code{scale_free}), a stochastic block
#' model (\code{stochastic_block_model}), or Watts–Strogatz small-world
#' (\code{small_world}) \insertCite{barabasiEmergenceScalingRandom1999}{DeCovarT}.
#'
#' @param n_genes Integer \eqn{G}, number of nodes (genes).
#' @param graph_model One of \code{"scale_free"},
#'   \code{"stochastic_block_model"}, \code{"small_world"}.
#' @param graph_params Named list of generator parameters:
#'   \describe{
#'     \item{\code{scale_free}}{\code{power}, \code{edges_per_node}
#'       (\code{m} in \code{igraph::sample_pa()})}
#'     \item{\code{stochastic_block_model}}{\code{block_prob},
#'       \code{p_within}, \code{p_between}}
#'     \item{\code{small_world}}{\code{nei}, \code{p} (rewiring
#'       probability) for \code{igraph::sample_smallworld()}}
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

  graph <- switch(
    graph_model,

    scale_free = {
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
      igraph::sample_smallworld(
        dim = 1L,
        size = n_genes,
        nei = nei,
        p = p_rewire
      )
    }
  )

  adjacency <- as.matrix(
    igraph::as_adjacency_matrix(igraph::as_undirected(graph))
  )
  storage.mode(adjacency) <- "integer"
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
#' from a random-graph model (\pkg{igraph}: scale-free, stochastic
#' block, or Watts–Strogatz small-world); i.i.d. signed weights
#' with inhibitory fraction \code{prop_inhibitory} form
#' \eqn{\boldsymbol{W}}; the precision is completed by a spectral shift;
#' each cell type shares
#' \eqn{\boldsymbol{\Sigma}_j=\boldsymbol{\Omega}^{-1}}.
#'
#' @param n_genes Integer; number of genes \eqn{G}.
#' @param n_celltypes Integer; number of cell types \eqn{J} (default 2).
#' @param mean_scale Positive scalar \eqn{s} for centroid norms
#'   (default \code{10}, as in the nine-scenario grid).
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
