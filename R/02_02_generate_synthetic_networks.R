# ---- simulate_synthetic_GRNs.R --------------------------------------------
##
## Generative model for Gaussian-based gene regulatory networks with
## AutoGeneS-inspired mean profiles and graph-constrained precision
## matrices.

#' Score mean profiles with AutoGeneS-style objectives
#'
#' Computes the two objectives used by AutoGeneS for signature quality.
#' For columns \eqn{\boldsymbol{\mu}_{\cdot j},\boldsymbol{\mu}_{\cdot k}}
#' of \eqn{\boldsymbol{\mu}}, the **cosine** (angle similarity, akin to
#' Pearson correlation of centred, unit-norm vectors but here without
#' mean-centring—only \eqn{\ell_2} normalisation) is
#' \deqn{
#'   \cos(\boldsymbol{\mu}_{\cdot j},\boldsymbol{\mu}_{\cdot k})
#'   =
#'   \frac{
#'     \boldsymbol{\mu}_{\cdot j}^{\mathsf{T}}
#'     \boldsymbol{\mu}_{\cdot k}
#'   }{
#'     \|\boldsymbol{\mu}_{\cdot j}\|_2\,
#'     \|\boldsymbol{\mu}_{\cdot k}\|_2
#'   },
#' }
#' and the **Euclidean** separation is
#' \deqn{
#'   \|\boldsymbol{\mu}_{\cdot j}-\boldsymbol{\mu}_{\cdot k}\|_2
#'   =
#'   \sqrt{
#'     \|\boldsymbol{\mu}_{\cdot j}\|_2^{2}
#'     +\|\boldsymbol{\mu}_{\cdot k}\|_2^{2}
#'     -2\,
#'     \boldsymbol{\mu}_{\cdot j}^{\mathsf{T}}
#'     \boldsymbol{\mu}_{\cdot k}
#'   }.
#' }
#' The returned scores are the mean absolute pairwise cosine (to minimise)
#' and the sum of pairwise Euclidean distances (to maximise).
#'
#' @param mean_signature_matrix Numeric matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (columns = cell types).
#'
#' @return A named list with \code{mean_abs_cosine} and
#'   \code{sum_euclidean_distance}.
#'
#' @keywords internal
#' @examples
#' compute_mean_profile_objectives(matrix(c(20, 22, 22, 20), 2))
#' @export
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
      # Cosine = <mu_j, mu_k> / (||mu_j|| ||mu_k||); 0 if a column is null.
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
#' The private vectors \eqn{\boldsymbol{v}_{j}} are
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
#' @details
#' **Private marker blocks.** Genes are partitioned into \eqn{J} nearly equal
#' contiguous blocks. Type \eqn{j}'s private direction \eqn{\boldsymbol{v}_{j}}
#' is the indicator of its block, then \eqn{\ell_2}-normalised. Distinct blocks
#' are orthogonal, so type-specific signal does not leak across columns before
#' the shared component is added.
#'
#' **Shared–private blend.** With unit shared direction
#' \eqn{\boldsymbol{u}=G^{-1/2}\mathbf{1}}, each column is
#' \eqn{\sqrt{\rho}\,\boldsymbol{u}+\sqrt{1-\rho}\,\boldsymbol{v}_{j}},
#' re-normalised, then scaled by \eqn{s}. Thus \eqn{\rho} dials collinearity
#' while \eqn{s} dials Euclidean separation without changing angles.
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
#' @examples
#' generate_mean_signature_matrix(
#'   n_genes = 6L,
#'   n_celltypes = 2L,
#'   target_cosine = 0.5
#' )
#' @export
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

  # Partition genes into J contiguous blocks; each block yields an orthogonal
  # private unit vector v_j (type-specific markers).
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
  # Blend shared u and private v_j, then unit-normalise and scale by s.
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
#' Draws undirected skeletons with \pkg{igraph}: Erdős–Rényi
#' (\code{erdos_renyi}), hub / star modules (\code{hub}), Barabási–Albert
#' preferential attachment (\code{scale_free}), a stochastic block model
#' (\code{stochastic_block_model}), or Watts–Strogatz small-world
#' (\code{small_world})
#' \insertCite{barabasiEmergenceScalingRandom1999}{DeCovarT}.
#'
#' @param n_genes Integer \eqn{G}, number of nodes (genes).
#' @param graph_model One of \code{"erdos_renyi"}, \code{"hub"},
#'   \code{"scale_free"}, \code{"stochastic_block_model"},
#'   \code{"small_world"} (also \code{"star"}, mapped to \code{"hub"}).
#'   Matching is case-insensitive.
#'
#' @srrstats {G2.3} Restricted character input (`graph_model`).
#' @srrstats {G2.3a} Validated via `.match_arg_ci()` (a `match.arg()`
#'   equivalent).
#' @srrstats {G2.3b} Matching is case-insensitive (`tolower()`).
#' @param graph_params Named list of generator parameters:
#'   \describe{
#'     \item{\code{erdos_renyi}}{\code{p} edge probability for
#'       \code{igraph::sample_gnp()} (default targets average degree
#'       about 2)}
#'     \item{\code{hub}}{\code{n_hubs}: partition nodes into that many
#'       groups and connect each group's hub to all other members (star
#'       modules)}
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
#' @examples
#' generate_random_network_skeleton(8L, "erdos_renyi")
#' @export
generate_random_network_skeleton <- function(
  n_genes,
  graph_model,
  graph_params = list()
) {
  n_genes <- as.integer(n_genes)
  graph_model <- .match_arg_ci(
    graph_model,
    c(
      "erdos_renyi",
      "hub",
      "star",
      "scale_free",
      "stochastic_block_model",
      "small_world"
    )
  )
  if (identical(graph_model, "star")) {
    graph_model <- "hub"
  }

  graph <- switch(
    graph_model,

    erdos_renyi = {
      p_edge <- if (is.null(graph_params$p)) {
        # Expected average degree ≈ 2 (sparse ER baseline).
        min(1, 2 / max(n_genes - 1L, 1L))
      } else {
        graph_params$p
      }
      igraph::sample_gnp(n_genes, p = p_edge, directed = FALSE)
    },

    hub = {
      n_hubs <- if (is.null(graph_params$n_hubs)) {
        max(1L, as.integer(round(sqrt(n_genes))))
      } else {
        as.integer(graph_params$n_hubs)
      }
      n_hubs <- max(1L, min(n_hubs, n_genes))
      # Equal-ish groups; first node of each group is the hub (star).
      block_sizes <- rep(n_genes %/% n_hubs, n_hubs)
      rem_hubs <- n_genes %% n_hubs
      if (rem_hubs > 0L) {
        block_sizes[seq_len(rem_hubs)] <-
          block_sizes[seq_len(rem_hubs)] + 1L
      }
      g <- igraph::make_empty_graph(n = n_genes, directed = FALSE)
      start <- 1L
      for (h in seq_len(n_hubs)) {
        members <- start:(start + block_sizes[[h]] - 1L)
        hub_node <- members[[1L]]
        leaves <- members[-1L]
        if (length(leaves) > 0L) {
          g <- igraph::add_edges(
            g,
            as.vector(rbind(hub_node, leaves))
          )
        }
        start <- start + block_sizes[[h]]
      }
      g
    },

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
#' @examples
#' A <- generate_random_network_skeleton(6L, "hub")
#' assign_iid_signed_weights(A, prop_inhibitory = 0.5)
#' @export
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
#' @examples
#' A <- generate_random_network_skeleton(6L, "hub")
#' W <- assign_iid_signed_weights(A)
#' build_normalised_precision(W, precision_shift = 0.2)
#' @export
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


#' Build cell-type-specific covariances from precision slices
#'
#' For each cell type \eqn{j}, inverts
#' \eqn{\boldsymbol{\Omega}_j} to
#' \eqn{\boldsymbol{\Sigma}_j=\boldsymbol{\Omega}_j^{-1}} and stacks the
#' result as a \eqn{G\times G\times J} array. Precisions need not be
#' shared across cell types.
#'
#' @param precision_array Symmetric positive-definite array
#'   \eqn{(\boldsymbol{\Omega}_j)_{j}\in\mathcal{M}_{G\times G\times J}}.
#'
#' @return Numeric array of dimension \eqn{G\times G\times J}.
#'
#' @keywords internal
#' @examples
#' A <- generate_random_network_skeleton(6L, "hub")
#' W <- assign_iid_signed_weights(A)
#' Omega <- build_normalised_precision(W, precision_shift = 0.2)
#' arr <- array(Omega, dim = c(6, 6, 1))
#' dim(build_covariance_array_from_precision(arr))
#' @export
build_covariance_array_from_precision <- function(precision_array) {
  if (
    is.null(dim(precision_array)) ||
      length(dim(precision_array)) != 3L
  ) {
    stop("`precision_array` must be a G x G x J array.", call. = FALSE)
  }
  n_genes <- dim(precision_array)[[1L]]
  n_celltypes <- dim(precision_array)[[3L]]
  if (dim(precision_array)[[2L]] != n_genes) {
    stop("`precision_array` dims must be G x G x J.", call. = FALSE)
  }

  gene_names <- dimnames(precision_array)[[1L]]
  celltype_names <- dimnames(precision_array)[[3L]]
  covariance_array <- array(
    NA_real_,
    dim = c(n_genes, n_genes, n_celltypes),
    dimnames = list(gene_names, gene_names, celltype_names)
  )
  for (j in seq_len(n_celltypes)) {
    sigma_j <- solve(precision_array[,, j])
    covariance_array[,, j] <- (sigma_j + t(sigma_j)) / 2
  }
  covariance_array
}


#' Simulate GRN first- and second-order moments
#'
#' @description
#' Builds a mean matrix \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}
#' and **cell-type-specific** second-order moments
#' \eqn{(\boldsymbol{\Omega}_j,\boldsymbol{\Sigma}_j=\boldsymbol{\Omega}_j^{-1})_{j=1}^{J}}
#' under a graph-constrained precision model. Means follow the
#' AutoGeneS-inspired cosine construction of
#' `generate_mean_signature_matrix()` with target pairwise cosine
#' \eqn{\rho}. For each cell type, an adjacency is drawn from a
#' random-graph model (or supplied), i.i.d. signed weights with
#' inhibitory fraction \code{prop_inhibitory} form \eqn{\boldsymbol{W}_j},
#' and the precision is completed by a spectral shift. Distinct cell
#' types receive **independent** precision draws by default (biology
#' rarely shares one network across types); pass a length-\eqn{J}
#' \code{graph_model} / \code{graph_params} or a pre-built
#' \code{adjacency} list / array for hybrid designs.
#'
#' @inheritParams generate_mean_signature_matrix
#' @param precision_shift Diagonal cushion \eqn{u} for the spectral
#'   shift (scalar, or length \eqn{J}).
#' @param precision_scale Positive magnitude \eqn{v} of signed
#'   off-diagonal precision weights (scalar, or length \eqn{J}).
#' @param prop_inhibitory Numeric in \eqn{[0,1]}; fraction of edges with
#'   positive precision weight (inhibitory partial correlation). Default
#'   \code{0.5} balances inhibitory and activatory edges (scalar, or
#'   length \eqn{J}).
#' @param graph_model One of \code{"erdos_renyi"}, \code{"hub"},
#'   \code{"scale_free"}, \code{"stochastic_block_model"},
#'   \code{"small_world"}, or a character vector of length \eqn{J}.
#' @param graph_params Named list of generator parameters (shared), or a
#'   list of length \eqn{J} of such named lists (see
#'   \code{generate_random_network_skeleton()}).
#' @param adjacency Optional pre-built undirected adjacencies: a list of
#'   \eqn{J} \eqn{G\times G} matrices, or a \eqn{G\times G\times J}
#'   array. When supplied, \code{graph_model} / \code{graph_params} are
#'   ignored for skeleton generation.
#'
#' @return Named list with:
#' * `mean_profiles`: matrix \eqn{\boldsymbol{\mu}};
#' * `covariance_matrices`: array
#'   \eqn{(\boldsymbol{\Sigma}_j)_{j}\in\mathcal{M}_{G\times G\times J}};
#' * `precision_matrices`: array
#'   \eqn{(\boldsymbol{\Omega}_j)_{j}\in\mathcal{M}_{G\times G\times J}};
#' * `graph_structure`: `adjacency_matrices`, `weighted_adjacencies`,
#'   and `normalised_precision` (all \eqn{G\times G\times J});
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
    "erdos_renyi",
    "hub",
    "star",
    "scale_free",
    "stochastic_block_model",
    "small_world"
  ),
  graph_params = list(),
  adjacency = NULL
) {
  n_genes <- as.integer(n_genes)
  n_celltypes <- as.integer(n_celltypes)

  if (n_genes < n_celltypes) {
    stop("`n_genes` must be at least `n_celltypes`.")
  }

  .recycle_positive <- function(x, nm) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, n_celltypes)
    }
    if (length(x) != n_celltypes || !all(is.finite(x)) || any(x <= 0)) {
      stop(
        "`",
        nm,
        "` must be a positive scalar or length-J vector.",
        call. = FALSE
      )
    }
    x
  }
  .recycle_unit <- function(x, nm) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      x <- rep(x, n_celltypes)
    }
    if (
      length(x) != n_celltypes ||
        !all(is.finite(x)) ||
        any(x < 0) ||
        any(x > 1)
    ) {
      stop(
        "`",
        nm,
        "` must lie in [0, 1] (scalar or length J).",
        call. = FALSE
      )
    }
    x
  }

  precision_shift <- .recycle_positive(precision_shift, "precision_shift")
  precision_scale <- .recycle_positive(precision_scale, "precision_scale")
  prop_inhibitory <- .recycle_unit(prop_inhibitory, "prop_inhibitory")

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

  ## --- step 2: per-cell-type adjacencies -------------------------------
  adjacency_array <- array(
    0L,
    dim = c(n_genes, n_genes, n_celltypes),
    dimnames = list(gene_names, gene_names, celltype_names)
  )

  if (!is.null(adjacency)) {
    if (is.list(adjacency)) {
      if (length(adjacency) != n_celltypes) {
        stop("`adjacency` list must have length J.", call. = FALSE)
      }
      for (j in seq_len(n_celltypes)) {
        a_j <- as.matrix(adjacency[[j]])
        if (!identical(dim(a_j), c(n_genes, n_genes))) {
          stop("Each `adjacency[[j]]` must be G x G.", call. = FALSE)
        }
        storage.mode(a_j) <- "integer"
        diag(a_j) <- 0L
        adjacency_array[,, j] <- a_j
      }
    } else {
      if (
        is.null(dim(adjacency)) ||
          length(dim(adjacency)) != 3L ||
          !identical(dim(adjacency), c(n_genes, n_genes, n_celltypes))
      ) {
        stop("`adjacency` array must be G x G x J.", call. = FALSE)
      }
      adjacency_array[] <- as.integer(adjacency)
      for (j in seq_len(n_celltypes)) {
        diag(adjacency_array[,, j]) <- 0L
      }
    }
  } else {
    model_choices <- c(
      "erdos_renyi",
      "hub",
      "star",
      "scale_free",
      "stochastic_block_model",
      "small_world"
    )
    if (length(graph_model) == 1L) {
      graph_model <- rep(
        .match_arg_ci(graph_model, model_choices),
        n_celltypes
      )
    } else {
      if (length(graph_model) != n_celltypes) {
        stop("`graph_model` must be length 1 or J.", call. = FALSE)
      }
      graph_model <- .match_arg_ci(graph_model, model_choices)
    }
    graph_model[graph_model == "star"] <- "hub"

    params_per_type <- if (
      length(graph_params) == n_celltypes &&
        all(vapply(graph_params, is.list, logical(1))) &&
        is.null(names(graph_params))
    ) {
      graph_params
    } else {
      rep(list(graph_params), n_celltypes)
    }

    for (j in seq_len(n_celltypes)) {
      a_j <- generate_random_network_skeleton(
        n_genes = n_genes,
        graph_model = graph_model[[j]],
        graph_params = params_per_type[[j]]
      )
      adjacency_array[,, j] <- a_j
    }
  }

  ## --- step 3: i.i.d. signs then spectral SPD completion, per type -----
  weighted_array <- array(
    0,
    dim = c(n_genes, n_genes, n_celltypes),
    dimnames = list(gene_names, gene_names, celltype_names)
  )
  precision_array <- array(
    NA_real_,
    dim = c(n_genes, n_genes, n_celltypes),
    dimnames = list(gene_names, gene_names, celltype_names)
  )

  for (j in seq_len(n_celltypes)) {
    weighted_array[,, j] <- assign_iid_signed_weights(
      adjacency_matrix = adjacency_array[,, j],
      prop_inhibitory = prop_inhibitory[[j]],
      weight_magnitude = precision_scale[[j]]
    )
    omega_j <- build_normalised_precision(
      weighted_adjacency = weighted_array[,, j],
      precision_shift = precision_shift[[j]]
    )
    precision_array[,, j] <- omega_j
  }

  covariance_matrices <- build_covariance_array_from_precision(
    precision_array
  )

  list(
    mean_profiles = mean_profiles,
    covariance_matrices = covariance_matrices,
    precision_matrices = precision_array,
    graph_structure = list(
      adjacency_matrices = adjacency_array,
      weighted_adjacencies = weighted_array,
      normalised_precision = precision_array
    ),
    objectives = objectives
  )
}
