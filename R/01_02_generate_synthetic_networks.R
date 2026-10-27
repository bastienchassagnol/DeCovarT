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


#' Thin orthonormal gene frame \eqn{\mathbf{Q}} for the mean signature
#'
#' With `nonnegative = TRUE`, the columns are indicators of contiguous,
#' near-equal gene blocks scaled to unit norm, so \eqn{\mathbf{Q}\ge 0} and
#' the means have disjoint, marker-like supports. Otherwise the frame is
#' the thin QR factor of a Gaussian draw (when `seed` is given; the
#' previous random state is restored) or of a fixed identity-plus-constant
#' matrix, which spreads the signal over all genes.
#'
#' @noRd
.gene_frame <- function(n_genes, n_celltypes, nonnegative, seed) {
  if (isTRUE(nonnegative)) {
    sizes <- rep(n_genes %/% n_celltypes, n_celltypes)
    extra <- n_genes %% n_celltypes
    sizes[seq_len(extra)] <- sizes[seq_len(extra)] + 1L
    q <- matrix(0, n_genes, n_celltypes)
    q[cbind(seq_len(n_genes), rep(seq_len(n_celltypes), sizes))] <-
      rep(1 / sqrt(sizes), sizes)
    return(q)
  }
  z <- if (!is.null(seed)) {
    .with_seed(seed, matrix(stats::rnorm(n_genes * n_celltypes), n_genes))
  } else {
    z <- matrix(1 / sqrt(n_genes), n_genes, n_celltypes)
    z[seq_len(n_celltypes), ] <- diag(n_celltypes)
    z
  }
  qr.Q(qr(z), complete = FALSE)[, seq_len(n_celltypes), drop = FALSE]
}

#' Generate mean profiles from a target Gram matrix
#'
#' Builds
#' \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} whose *unit* columns
#' realise a prescribed Gram matrix \eqn{R} (pairwise cosines) via the
#' symmetric square root,
#' \deqn{
#'   \boldsymbol{\mu}
#'   =
#'   s\,
#'   \boldsymbol{Q}
#'   R^{1/2},
#'   \qquad
#'   \boldsymbol{Q}^{\mathsf{T}}\boldsymbol{Q}=\boldsymbol{I}_{J}.
#' }
#' Then \eqn{\boldsymbol{\mu}^{\mathsf{T}}\boldsymbol{\mu}=s^{2}R} whenever
#' \eqn{R} has unit diagonal. The default \eqn{R} is the equicorrelation
#' (equicorrelation) matrix
#' \eqn{(1-\rho)I+\rho\mathbf{1}\mathbf{1}^{\mathsf{T}}}, positive
#' semidefinite for \eqn{\rho\in[-1/(J-1),1]}. Supply `target_gram` to set unequal
#' pairwise cosines (for example two close pairs and a distant background),
#' provided \(R\) is symmetric nonnegative-definite.
#'
#' Realisations are unique up to an orthogonal transformation of the gene
#' space: replacing \eqn{\boldsymbol{Q}} by \eqn{\boldsymbol{Q}U} with
#' \eqn{U^{\mathsf{T}}U=\boldsymbol{I}_{J}} leaves the Gram unchanged. By
#' default \eqn{\boldsymbol{Q}} is a deterministic thin QR frame; pass `seed`
#' for a Haar-like Gaussian frame (the global random state is restored). Cholesky \eqn{R=LL^{\mathsf{T}}} is an
#' alternative square root used in Monte Carlo simulation; see the
#' how-to-build-synthetic-scenarios vignette appendix.
#'
#' @param n_genes Integer \eqn{G}; must be at least \code{n_celltypes}.
#' @param n_celltypes Integer \eqn{J\ge 2}.
#' @param mean_scale Positive scalar \eqn{s} (centroid norms). Default
#'   \code{10}. Hold fixed when studying cosine / collinearity alone.
#' @param target_cosine Numeric pairwise cosine \eqn{\rho} for the
#'   equicorrelation Gram. Ignored when `target_gram` is supplied.
#'   Must lie in \eqn{[-1/(J-1),1]}.
#' @param target_gram Optional symmetric \eqn{J\times J} cosine matrix
#'   (unit diagonal). Overrides `target_cosine`.
#' @param seed Optional integer. When supplied, the orthonormal frame
#'   \eqn{\boldsymbol{Q}} is drawn from a Gaussian QR; otherwise it is
#'   deterministic. Ignored when `nonnegative = TRUE`.
#' @param nonnegative If `TRUE`, use a disjoint-support nonnegative
#'   frame so that \eqn{\boldsymbol{\mu}\ge 0} whenever \eqn{R^{1/2}\ge 0}
#'   (required by [deconvolute_ratios()]). The Gram is unchanged.
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
#' k <- matrix(c(1, 0.98, 0.2, 0.98, 1, 0.2, 0.2, 0.2, 1), 3, 3)
#' generate_mean_signature_matrix(
#'   n_genes = 8L,
#'   n_celltypes = 3L,
#'   target_gram = k
#' )
#' @export
generate_mean_signature_matrix <- function(
  n_genes,
  n_celltypes,
  mean_scale = 10,
  target_cosine = 0,
  target_gram = NULL,
  seed = NULL,
  nonnegative = FALSE,
  gene_names = NULL,
  celltype_names = NULL
) {
  n_genes <- as.integer(n_genes)
  n_celltypes <- as.integer(n_celltypes)
  if (length(n_genes) != 1L || is.na(n_genes) || n_genes < n_celltypes) {
    stop("`n_genes` must be at least `n_celltypes`.", call. = FALSE)
  }
  if (
    length(n_celltypes) != 1L ||
      is.na(n_celltypes) ||
      n_celltypes < 2L
  ) {
    stop("`n_celltypes` must be at least 2.", call. = FALSE)
  }
  if (mean_scale <= 0) {
    stop("`mean_scale` must be positive.", call. = FALSE)
  }

  if (is.null(target_gram)) {
    rho <- as.numeric(target_cosine)
    rho_min <- -1 / (n_celltypes - 1)
    if (
      length(rho) != 1L ||
        is.na(rho) ||
        rho < rho_min - 1e-10 ||
        rho > 1 + 1e-10
    ) {
      stop(
        "`target_cosine` must be one number in [",
        format(rho_min, digits = 4),
        ", 1] for J = ",
        n_celltypes,
        ".",
        call. = FALSE
      )
    }
    gram <- (1 - rho) * diag(n_celltypes) + rho
  } else {
    gram <- as.matrix(target_gram)
    if (
      nrow(gram) != n_celltypes ||
        ncol(gram) != n_celltypes
    ) {
      stop(
        "`target_gram` must be a square J x J matrix.",
        call. = FALSE
      )
    }
    if (max(abs(gram - t(gram))) > 1e-8) {
      stop("`target_gram` must be symmetric.", call. = FALSE)
    }
    gram <- .symmetrise(gram)
    evals <- eigen(gram, symmetric = TRUE, only.values = TRUE)$values
    if (min(evals) < -1e-8) {
      stop(
        "`target_gram` must be nonnegative-definite.",
        call. = FALSE
      )
    }
  }

  q <- .gene_frame(n_genes, n_celltypes, nonnegative, seed)
  # Symmetric square root of the Gram, so that crossprod(mu) = s^2 R.
  eig <- eigen(gram, symmetric = TRUE)
  root <- eig$vectors %*%
    diag(sqrt(pmax(eig$values, 0)), n_celltypes) %*%
    t(eig$vectors)
  mean_signature_matrix <- mean_scale * q %*% root
  if (isTRUE(nonnegative)) {
    if (any(mean_signature_matrix < -1e-10)) {
      stop(
        "`nonnegative = TRUE` needs a Gram whose symmetric square root ",
        "is nonnegative (for example all pairwise cosines in [0, 1]).",
        call. = FALSE
      )
    }
    mean_signature_matrix[mean_signature_matrix < 0] <- 0
  }

  if (is.null(gene_names)) {
    gene_names <- paste0("gene_", seq_len(n_genes))
  }
  if (is.null(celltype_names)) {
    celltype_names <- paste0("celltype_", seq_len(n_celltypes))
  }
  dimnames(mean_signature_matrix) <- list(gene_names, celltype_names)
  mean_signature_matrix
}


#' Sample a symmetric adjacency from a random-graph family
#'
#' Draws undirected skeletons with \pkg{igraph}: Erdős–Rényi
#' (\code{erdos_renyi}), hub / star modules (\code{hub}), Barabási–Albert
#' preferential attachment (\code{scale_free}), a stochastic block model
#' (\code{stochastic_block_model}), or Watts–Strogatz small-world
#' (\code{small_world})
#' \insertCite{barabasiEmergenceScalingRandom1999,hollandStochasticBlockmodelsFirst1983,wattsCollectiveDynamicsSmallworld1998}{DeCovarT}.
#'
#' @param n_genes Integer \eqn{G}, number of nodes (genes).
#' @param graph_model One of \code{"erdos_renyi"}, \code{"hub"},
#'   \code{"scale_free"}, \code{"stochastic_block_model"},
#'   \code{"small_world"} (also \code{"star"}, mapped to \code{"hub"}).
#'   Matching is case-insensitive.
#'
#' @srrstats {G2.3} Restricted character input (`graph_model`).
#' @srrstats {G2.3a} Validated via `.match_arg_case_insensitive()` (a
#'   `match.arg()` equivalent).
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
#'       \code{p_within}, \code{p_between}; or \code{block_sizes}, an
#'       integer vector summing to \code{n_genes} that fixes the blocks
#'       (contiguous genes) instead of drawing them from
#'       \code{block_prob}}
#'     \item{\code{small_world}}{\code{nei}, \code{p} (rewiring
#'       probability) for \code{igraph::sample_smallworld()}}
#'   }
#'
#' @return Symmetric integer matrix \eqn{G\times G} with zero diagonal.
#'   For a stochastic block model, attribute \code{"block"} gives the
#'   block of every gene.
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
  .check_suggested_package("igraph", "generate_random_network_skeleton")
  n_genes <- as.integer(n_genes)
  graph_model <- .match_arg_case_insensitive(
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

      if (is.null(graph_params$block_sizes)) {
        block_sizes <- as.integer(
          c(stats::rmultinom(1L, n_genes, block_prob))
        )
        block_sizes <- pmax(block_sizes, 1L)
        size_diff <- n_genes - sum(block_sizes)
        if (size_diff != 0L) {
          idx <- which.max(block_sizes)
          block_sizes[idx] <- block_sizes[idx] + size_diff
        }
      } else {
        block_sizes <- as.integer(graph_params$block_sizes)
        if (sum(block_sizes) != n_genes || any(block_sizes < 1L)) {
          stop(
            "`block_sizes` must be positive and sum to `n_genes`.",
            call. = FALSE
          )
        }
        n_blocks <- length(block_sizes)
        pref_matrix <- matrix(p_between, n_blocks, n_blocks)
        diag(pref_matrix) <- p_within
      }
      block <- rep(seq_along(block_sizes), block_sizes)

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
  if (identical(graph_model, "stochastic_block_model")) {
    attr(adjacency, "block") <- block
  }
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
  if (prop_inhibitory < 0 || prop_inhibitory > 1) {
    stop("`prop_inhibitory` must lie in [0, 1].")
  }
  if (weight_magnitude <= 0) {
    stop("`weight_magnitude` must be positive.")
  }

  W <- matrix(0, nrow(A), ncol(A))
  # Non-null off-diagonal support of A (upper triangle). If every off-diagonal
  # entry is zero there is no edge whose sign can change, so return W = 0.
  upper <- which(upper.tri(A) & A != 0, arr.ind = TRUE)
  if (nrow(upper) == 0L) {
    dimnames(W) <- dimnames(A)
    return(W)
  }

  n_edge <- nrow(upper)
  n_inhib <- as.integer(round(prop_inhibitory * n_edge))
  signs <- rep(-1, n_edge) # activatory precision weight by default
  # Sample which edges receive a negative partial correlation (inhibitory
  # precision weight): uniform without replacement over the edge set.
  if (n_inhib > 0L) {
    inhib_idx <- sample.int(n_edge, n_inhib)
    signs[inhib_idx] <- 1
  }
  values <- signs * weight_magnitude
  # Symmetrise: undirected Gaussian graphical models require a symmetric
  # precision (or covariance) matrix, so Omega_jk = Omega_kj.
  W[upper] <- values
  W[cbind(upper[, 2L], upper[, 1L])] <- values
  dimnames(W) <- dimnames(A)
  W
}


#' Uniform spectral shift that preserves off-diagonal support
#'
#' Completes a symmetric matrix by
#' \eqn{A\leftarrow A+(\lvert\min(0,\lambda_{\min})\rvert+u)\mathbf{I}}.
#' Off-diagonals (graph support and signs) are unchanged.
#'
#' @keywords internal
#' @noRd
.uniform_spectral_shift <- function(mat, cushion) {
  mat <- .symmetrise(mat)
  cushion <- as.numeric(cushion)
  if (length(cushion) != 1L || !is.finite(cushion) || cushion <= 0) {
    stop("`cushion` must be a positive scalar.", call. = FALSE)
  }
  evals <- eigen(mat, symmetric = TRUE, only.values = TRUE)$values
  extra <- max(0, -min(evals)) + cushion
  diag(mat) <- diag(mat) + extra
  mat
}

#' Repeat a uniform spectral shift until Cholesky succeeds
#'
#' @keywords internal
#' @noRd
.ensure_spd <- function(
  mat,
  cushion = sqrt(.Machine$double.eps),
  max_tries = 8L
) {
  mat <- .symmetrise(mat)
  for (k in seq_len(as.integer(max_tries))) {
    if (!inherits(try(chol(mat), silent = TRUE), "try-error")) {
      return(mat)
    }
    mat <- .uniform_spectral_shift(
      mat,
      cushion = cushion * (2^(k - 1L))
    )
  }
  stop(
    "Could not complete a symmetric matrix to strictly positive ",
    "definite form by a uniform spectral shift.",
    call. = FALSE
  )
}

#' Invert an SPD precision to a Cholesky-able covariance
#'
#' \eqn{\boldsymbol{\Sigma}=\boldsymbol{\Omega}^{-1}} from the Cholesky
#' factor of \eqn{\boldsymbol{\Omega}} ([chol2inv()]). For a symmetric
#' positive-definite precision this is exactly symmetric, at least as
#' accurate as a QR inverse and three to five times faster, and the
#' factorisation itself certifies positive definiteness. If the factor of
#' \eqn{\boldsymbol{\Sigma}} fails
#' (an ill-conditioned hub graph), the support-preserving spectral shift
#' of \eqn{\boldsymbol{\Omega}} is doubled until both factors exist.
#'
#' @return A list with the (possibly shifted) `precision` and its
#'   `covariance`.
#'
#' @keywords internal
#' @noRd
.invert_precision <- function(omega, cushion) {
  omega <- .ensure_spd(omega, cushion = cushion)
  extra <- cushion
  for (k in seq_len(8L)) {
    sigma <- chol2inv(chol(omega))
    if (!inherits(try(chol(sigma), silent = TRUE), "try-error")) {
      dimnames(sigma) <- dimnames(omega)
      return(list(precision = omega, covariance = sigma))
    }
    extra <- extra * 2
    omega <- .uniform_spectral_shift(omega, cushion = extra)
  }
  stop(
    "Could not invert a precision to a strictly positive-definite ",
    "covariance by a uniform spectral shift.",
    call. = FALSE
  )
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
#'   \bigl(\lvert\min(0,\lambda_{\min}(\boldsymbol{W}))\rvert + u\bigr)
#'   \mathbf{I}.
#' }
#' The result is checked with [base::chol()]; if a factor still fails
#' (ill-conditioned graphs), the same loading is increased until
#' \eqn{\boldsymbol{\Omega}\succ 0} strictly.
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
  omega <- .symmetrise(weighted_adjacency)
  diag(omega) <- 0
  omega <- .uniform_spectral_shift(omega, cushion = precision_shift)
  .ensure_spd(omega, cushion = precision_shift)
}


#' Build cell-type-specific covariances from precision slices
#'
#' For each cell type \eqn{j}, inverts
#' \eqn{\boldsymbol{\Omega}_j} to
#' \eqn{\boldsymbol{\Sigma}_j=\boldsymbol{\Omega}_j^{-1}} (via
#' [chol2inv()]) and stacks the result as a \eqn{G\times G\times J}
#' array. If a slice is not numerically SPD, or if the inverse fails
#' Cholesky, a further uniform spectral shift is applied to that
#' precision (off-diagonal support unchanged). Precisions need not be
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
    covariance_array[,, j] <- .invert_precision(
      precision_array[,, j],
      cushion = sqrt(.Machine$double.eps)
    )$covariance
  }
  covariance_array
}


#' Simulate GRN first- and second-order moments
#'
#' @description
#' Builds a mean matrix \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}}
#' and **cell-type-specific** second-order moments
#' \eqn{(\boldsymbol{\Omega}_j,
#'   \boldsymbol{\Sigma}_j=\boldsymbol{\Omega}_j^{-1})_{j=1}^{J}}
#' under a graph-constrained precision model. Means realise a target
#' Gram matrix via the symmetric square root of
#' `generate_mean_signature_matrix()` (equicorrelation \eqn{\rho} by
#' default, or a full `target_gram`). For each cell type, an adjacency is drawn from a
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
  target_gram = NULL,
  seed = NULL,
  nonnegative = FALSE,
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

  ## --- step 1: Gram-matrix mean profiles -------------------------------
  mean_profiles <- generate_mean_signature_matrix(
    n_genes = n_genes,
    n_celltypes = n_celltypes,
    mean_scale = mean_scale,
    target_cosine = target_cosine,
    target_gram = target_gram,
    seed = seed,
    nonnegative = nonnegative,
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
        .match_arg_case_insensitive(graph_model, model_choices),
        n_celltypes
      )
    } else {
      if (length(graph_model) != n_celltypes) {
        stop("`graph_model` must be length 1 or J.", call. = FALSE)
      }
      graph_model <- .match_arg_case_insensitive(graph_model, model_choices)
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
  covariance_array <- array(
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
    pair <- .invert_precision(
      build_normalised_precision(weighted_array[,, j], precision_shift[[j]]),
      cushion = precision_shift[[j]]
    )
    precision_array[,, j] <- pair$precision
    covariance_array[,, j] <- pair$covariance
  }

  list(
    mean_profiles = mean_profiles,
    covariance_matrices = covariance_array,
    precision_matrices = precision_array,
    graph_structure = list(
      adjacency_matrices = adjacency_array,
      weighted_adjacencies = weighted_array,
      normalised_precision = precision_array
    ),
    objectives = objectives
  )
}
