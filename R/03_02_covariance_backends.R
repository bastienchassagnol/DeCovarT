# Structure-aware factorisation backends for the bulk covariance
# Sigma(p) = sum_j p_j^2 Sigma_j. The dense Cholesky in
# .sigma_p_factorisation() is the universal default; the backends here
# exploit *declared* structure (block, band, sparse, diagonal-plus-low-rank)
# to obtain log-determinant and solves without a dense G x G inverse.
#
# Design follows the numerical-linear-algebra notes: expose OPERATORS
# (logdet, solve, quadform, trace of precision times a matrix) rather than a
# materialised precision, so a sparse or Woodbury solve is not silently
# undone by a dense chol2inv. Structure is supplied explicitly, never
# inferred from floating-point near-zeros.

#' Recommend a covariance backend from a network topology
#'
#' @description
#' Maps a `graph_model` name (see
#' [generate_random_network_skeleton()]) to the covariance-factorisation
#' `structure` most likely to help. This encodes the decision-tree workflow:
#' the dense Cholesky is always a safe default; the structured backends are
#' refinements when the covariance (not merely the precision) inherits
#' exploitable structure.
#'
#' @details
#' A sparse gene-network **precision** matrix
#' \eqn{\boldsymbol{\Omega}_j=\boldsymbol{\Sigma}_j^{-1}} does **not** imply a
#' sparse **covariance** \eqn{\boldsymbol{\Sigma}_j}; the inverse of a sparse
#' precision is generally dense. The mapping below therefore expresses
#' modelling intent, and the returned structure only accelerates the
#' likelihood when \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} itself (or a
#' cheap representation of it) has that structure. Disconnected gene modules
#' give an exactly block-diagonal covariance; shared regulatory programs give
#' a diagonal-plus-low-rank covariance; a small set of bridge / housekeeping
#' hubs gives an arrow pattern handled by the Schur complement.
#'
#' @param graph_model One of `"erdos_renyi"`, `"hub"`, `"star"`,
#'   `"scale_free"`, `"stochastic_block_model"`, `"small_world"`, `"band"`,
#'   `"ar"`. Matching is case-insensitive.
#'
#' @return A single string: one of `"dense"`, `"block"`, `"band"`,
#'   `"sparse"`, `"diag_lowrank"`.
#'
#' @keywords internal
#' @examples
#' covariance_structure_from_graph_model("stochastic_block_model")
#' covariance_structure_from_graph_model("erdos_renyi")
#' @seealso [new_decovart_covariance()]
#' @export
covariance_structure_from_graph_model <- function(graph_model) {
  graph_model <- .match_arg_case_insensitive(
    graph_model,
    c(
      "erdos_renyi",
      "hub",
      "star",
      "scale_free",
      "stochastic_block_model",
      "small_world",
      "band",
      "ar"
    )
  )
  switch(
    graph_model,
    stochastic_block_model = "block",
    hub = "diag_lowrank",
    star = "diag_lowrank",
    scale_free = "sparse",
    small_world = "sparse",
    erdos_renyi = "sparse",
    band = "band",
    ar = "band"
  )
}

#' Structure-aware covariance backend for \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'
#' @description
#' Wraps the cell-type covariance array together with a declared
#' `structure`, so that the bulk covariance
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^2\boldsymbol{\Sigma}_j}
#' can be factorised with the cheapest exact method. The object exposes
#' operators via [sigma_logdet()], [sigma_solve()], [sigma_quadform()] and
#' [sigma_trace_precision_times()]; the dense array is always accepted and
#' wrapped as `"dense"`.
#'
#' @details
#' Supported structures:
#' * `"dense"`: universal fallback; one Cholesky of the assembled matrix.
#' * `"block"`: exactly block-diagonal covariance (disconnected gene
#'   modules); one Cholesky per block, so the factorisation cost drops by
#'   roughly the squared number of equal blocks.
#' * `"band"`: covariance with bandwidth \eqn{b}
#'   (\eqn{\Sigma_{jk}=0} for \eqn{|j-k|>b}); a banded / sparse Cholesky.
#' * `"sparse"`: sparse covariance with a fixed nonzero pattern across
#'   \eqn{\boldsymbol{p}}; the fill-reducing ordering (symbolic
#'   factorisation) is computed once and reused for every numeric refactor.
#' * `"diag_lowrank"`: \eqn{\boldsymbol{\Sigma}_j=\mathbf{D}_j+
#'   \mathbf{U}\mathbf{C}_j\mathbf{U}^{\mathsf{T}}} with shared loadings
#'   \eqn{\mathbf{U}}; solves use the Woodbury identity and the
#'   log-determinant uses the matrix-determinant lemma, transferring the
#'   \eqn{O(G^3)} cost to the low rank \eqn{r\ll G}.
#'
#' Band and sparse factorisations use \pkg{Matrix} (CHOLMOD;
#' \insertCite{chenAlgorithm887CHOLMOD2008}{DeCovarT}).
#' Dense Cholesky complexity follows
#' \insertCite{golubMatrixComputations2013}{DeCovarT};
#' the low-rank path uses the Woodbury identity
#' \insertCite{hagerUpdatingInverseMatrix1989}{DeCovarT}.
#'
#' @param Sigma Numeric array of cell-type covariances in
#'   \eqn{\mathcal{M}_{G\times G\times J}}. For `"diag_lowrank"` it may be
#'   omitted and is then reconstructed from `diagonal`, `loadings`, `core`.
#' @param structure One of `"dense"`, `"block"`, `"band"`, `"sparse"`,
#'   `"diag_lowrank"`. Matching is case-insensitive.
#' @param blocks For `"block"`: an integer / factor vector of length
#'   \eqn{G} giving block membership, or a list of index vectors.
#' @param bandwidth For `"band"`: non-negative integer bandwidth \eqn{b}.
#' @param diagonal For `"diag_lowrank"`: numeric \eqn{G\times J} matrix whose
#'   column \eqn{j} is the diagonal of \eqn{\mathbf{D}_j}.
#' @param loadings For `"diag_lowrank"`: shared loadings
#'   \eqn{\mathbf{U}\in\mathcal{M}_{G\times r}}.
#' @param core For `"diag_lowrank"`: array \eqn{\mathcal{M}_{r\times r\times J}}
#'   of cores \eqn{\mathbf{C}_j}.
#' @param tol Relative tolerance for validating that off-structure entries
#'   vanish (block / band).
#'
#' @return An object of class `"decovart_covariance"`.
#'
#' @keywords internal
#' @examples
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' cov_dense <- new_decovart_covariance(Sigma, "dense")
#' sigma_logdet(cov_dense, c(0.6, 0.4))
#' @seealso [sigma_logdet()], [covariance_structure_from_graph_model()]
#' @references
#' \insertRef{chenAlgorithm887CHOLMOD2008}{DeCovarT}
#'
#' \insertRef{golubMatrixComputations2013}{DeCovarT}
#'
#' \insertRef{hagerUpdatingInverseMatrix1989}{DeCovarT}
#' @export
new_decovart_covariance <- function(
  Sigma = NULL,
  structure = c("dense", "block", "band", "sparse", "diag_lowrank"),
  blocks = NULL,
  bandwidth = NULL,
  diagonal = NULL,
  loadings = NULL,
  core = NULL,
  tol = 1e-8
) {
  structure <- .match_arg_case_insensitive(
    structure,
    c("dense", "block", "band", "sparse", "diag_lowrank")
  )

  if (identical(structure, "diag_lowrank") && is.null(Sigma)) {
    Sigma <- .assemble_diag_lowrank_array(diagonal, loadings, core)
  }
  .validate_sigma_array(Sigma)

  n_genes <- dim(Sigma)[[1L]]
  n_celltypes <- dim(Sigma)[[3L]]

  backend <- list(
    Sigma = Sigma,
    structure = structure,
    n_genes = n_genes,
    n_celltypes = n_celltypes,
    cache = new.env(parent = emptyenv())
  )

  if (identical(structure, "block")) {
    backend$blocks <- .validate_covariance_blocks(Sigma, blocks, tol)
  } else if (identical(structure, "band")) {
    backend$bandwidth <- .validate_covariance_bandwidth(Sigma, bandwidth, tol)
  } else if (identical(structure, "diag_lowrank")) {
    lr <- .validate_diag_lowrank(Sigma, diagonal, loadings, core, tol)
    backend$diagonal <- lr$diagonal
    backend$loadings <- lr$loadings
    backend$core <- lr$core
  }

  structure(backend, class = "decovart_covariance")
}

#' @rdname new_decovart_covariance
#' @param x A `decovart_covariance` object (print method).
#' @param ... Ignored.
#' @export
#' @method print decovart_covariance
print.decovart_covariance <- function(x, ...) {
  cat(
    "<decovart_covariance>",
    "\n  structure:  ",
    x$structure,
    "\n  genes:      ",
    x$n_genes,
    "\n  cell types: ",
    x$n_celltypes,
    "\n"
  )
  if (identical(x$structure, "block")) {
    cat("  blocks:     ", length(x$blocks), "\n")
  } else if (identical(x$structure, "band")) {
    cat("  bandwidth:  ", x$bandwidth, "\n")
  } else if (identical(x$structure, "diag_lowrank")) {
    cat("  rank:       ", ncol(x$loadings), "\n")
  }
  invisible(x)
}

# Operators ---------------------------------------------------------------

#' Log-determinant of \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'
#' @param backend A [new_decovart_covariance()] object.
#' @param p Numeric proportions vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#'
#' @return Scalar \eqn{\log\det\boldsymbol{\Sigma}(\boldsymbol{p})}.
#'
#' @keywords internal
#' @examples
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' sigma_logdet(new_decovart_covariance(Sigma, "dense"), c(0.6, 0.4))
#' @seealso [sigma_solve()], [sigma_quadform()]
#' @export
sigma_logdet <- function(backend, p) {
  .sigma_prepare(backend, p)$logdet
}

#' Precision solve \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{B}}
#'
#' @inheritParams sigma_logdet
#' @param b Numeric vector or matrix of right-hand sides.
#'
#' @return \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{b}}, shaped
#'   like `b`.
#'
#' @keywords internal
#' @examples
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' sigma_solve(new_decovart_covariance(Sigma, "dense"), c(0.6, 0.4), c(1, 2))
#' @seealso [sigma_logdet()]
#' @export
sigma_solve <- function(backend, p, b) {
  .sigma_prepare(backend, p)$solve(b)
}

#' Mahalanobis quadratic form
#' \eqn{\mathbf{r}^{\mathsf{T}}\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}\mathbf{r}}
#'
#' @inheritParams sigma_logdet
#' @param r Numeric residual vector.
#'
#' @return Scalar quadratic form.
#'
#' @keywords internal
#' @examples
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' sigma_quadform(new_decovart_covariance(Sigma, "dense"), c(0.6, 0.4), c(1, 2))
#' @seealso [sigma_logdet()]
#' @export
sigma_quadform <- function(backend, p, r) {
  .sigma_prepare(backend, p)$quadform(r)
}

#' Trace of the precision times a matrix,
#' \eqn{\operatorname{tr}(\boldsymbol{\Theta}(\boldsymbol{p})\mathbf{S})}
#'
#' @description
#' Computes \eqn{\operatorname{tr}(\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
#' \mathbf{S})} through solves only, never materialising the dense precision.
#' This is the trace term that appears in the analytic gradient and Hessian.
#'
#' @inheritParams sigma_logdet
#' @param S Numeric \eqn{G\times G} matrix (typically a cell-type covariance
#'   slice \eqn{\boldsymbol{\Sigma}_j}).
#'
#' @return Scalar trace.
#'
#' @keywords internal
#' @examples
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' cov <- new_decovart_covariance(Sigma, "dense")
#' sigma_trace_precision_times(cov, c(0.6, 0.4), Sigma[, , 1])
#' @seealso [sigma_solve()]
#' @export
sigma_trace_precision_times <- function(backend, p, S) {
  solved <- .sigma_prepare(backend, p)$solve(S)
  sum(diag(solved))
}

# Preparation and dispatch ------------------------------------------------

#' Prepare (and cache) a structured factorisation at one trial p
#' @return List with `logdet`, `solve`, `quadform`, `matrix`, `chol`.
#' @keywords internal
#' @noRd
.sigma_prepare <- function(backend, p) {
  if (!inherits(backend, "decovart_covariance")) {
    stop("`backend` must be a `decovart_covariance` object.", call. = FALSE)
  }
  p <- as.numeric(p)
  if (length(p) != backend$n_celltypes) {
    stop(
      "`p` must have length ",
      backend$n_celltypes,
      ".",
      call. = FALSE
    )
  }

  cache <- backend$cache
  if (!is.null(cache$p) && identical(cache$p, p)) {
    return(cache$prep)
  }

  prep <- switch(
    backend$structure,
    dense = .sigma_prepare_dense(backend, p),
    block = .sigma_prepare_block(backend, p),
    band = .sigma_prepare_sparse(backend, p),
    sparse = .sigma_prepare_sparse(backend, p),
    diag_lowrank = .sigma_prepare_diag_lowrank(backend, p)
  )

  cache$p <- p
  cache$prep <- prep
  prep
}

#' Dense Cholesky (universal default)
#' @keywords internal
#' @noRd
.sigma_prepare_dense <- function(backend, p) {
  sigma_p <- .compute_global_variance(p, backend$Sigma)
  chol_factor <- chol(sigma_p)
  list(
    logdet = 2 * sum(log(diag(chol_factor))),
    solve = function(b) .solve_chol(chol_factor, b),
    quadform = function(r) {
      z <- backsolve(chol_factor, r, transpose = TRUE)
      sum(z * z)
    },
    matrix = sigma_p,
    chol = chol_factor
  )
}

#' Block Cholesky for an exactly block-diagonal covariance
#' @keywords internal
#' @noRd
.sigma_prepare_block <- function(backend, p) {
  sigma_p <- .compute_global_variance(p, backend$Sigma)
  blocks <- backend$blocks
  chols <- lapply(blocks, function(idx) chol(sigma_p[idx, idx, drop = FALSE]))
  logdet <- sum(vapply(
    chols,
    function(r) 2 * sum(log(diag(r))),
    numeric(1L)
  ))

  solve_fun <- function(b) {
    b_mat <- as.matrix(b)
    out <- matrix(0, nrow = nrow(b_mat), ncol = ncol(b_mat))
    for (k in seq_along(blocks)) {
      idx <- blocks[[k]]
      out[idx, ] <- .solve_chol(chols[[k]], b_mat[idx, , drop = FALSE])
    }
    if (is.null(dim(b))) out[, 1L] else out
  }

  list(
    logdet = logdet,
    solve = solve_fun,
    quadform = function(r) {
      total <- 0
      for (k in seq_along(blocks)) {
        idx <- blocks[[k]]
        z <- backsolve(chols[[k]], r[idx], transpose = TRUE)
        total <- total + sum(z * z)
      }
      total
    },
    matrix = sigma_p,
    chol = NULL
  )
}

#' Sparse / banded Cholesky with cached symbolic factorisation
#' @keywords internal
#' @noRd
.sigma_prepare_sparse <- function(backend, p) {
  sigma_p <- .compute_global_variance(p, backend$Sigma)
  sparse_sigma <- Matrix::forceSymmetric(
    Matrix::Matrix(sigma_p, sparse = TRUE)
  )

  cache <- backend$cache
  chol_factor <- if (!is.null(cache$symbolic)) {
    # Reuse the fill-reducing ordering computed once (symbolic step).
    tryCatch(
      Matrix::update(cache$symbolic, sparse_sigma),
      error = function(e) {
        Matrix::Cholesky(sparse_sigma, perm = TRUE, LDL = FALSE)
      }
    )
  } else {
    Matrix::Cholesky(sparse_sigma, perm = TRUE, LDL = FALSE)
  }
  cache$symbolic <- chol_factor

  logdet <- as.numeric(
    Matrix::determinant(sparse_sigma, logarithm = TRUE)$modulus
  )

  solve_fun <- function(b) {
    solved <- Matrix::solve(chol_factor, b, system = "A")
    solved <- as.matrix(solved)
    if (is.null(dim(b))) solved[, 1L] else solved
  }

  list(
    logdet = logdet,
    solve = solve_fun,
    quadform = function(r) sum(r * solve_fun(r)),
    matrix = sigma_p,
    chol = NULL
  )
}

#' Diagonal-plus-low-rank (Woodbury / matrix-determinant lemma)
#' @keywords internal
#' @noRd
.sigma_prepare_diag_lowrank <- function(backend, p) {
  p2 <- p^2
  d_vec <- as.numeric(backend$diagonal %*% p2)
  u_mat <- backend$loadings
  core_p <- .weighted_core_sum(backend$core, p2)

  d_inv <- 1 / d_vec
  # M = C(p)^{-1} + U^T D^{-1} U  (r x r, r << G)
  m_mat <- solve(core_p) + crossprod(u_mat, d_inv * u_mat)
  m_inv <- solve(m_mat)

  logdet <- sum(log(d_vec)) +
    as.numeric(determinant(core_p, logarithm = TRUE)$modulus) +
    as.numeric(determinant(m_mat, logarithm = TRUE)$modulus)

  solve_fun <- function(b) {
    b_mat <- as.matrix(b)
    dinv_b <- d_inv * b_mat
    correction <- d_inv * (u_mat %*% (m_inv %*% crossprod(u_mat, dinv_b)))
    out <- dinv_b - correction
    if (is.null(dim(b))) out[, 1L] else out
  }

  list(
    logdet = logdet,
    solve = solve_fun,
    quadform = function(r) sum(r * solve_fun(r)),
    matrix = NULL,
    chol = NULL
  )
}

# Helpers -----------------------------------------------------------------

#' Solve S X = B given the upper Cholesky factor R (S = R^T R)
#' @keywords internal
#' @noRd
.solve_chol <- function(chol_factor, b) {
  backsolve(
    chol_factor,
    backsolve(chol_factor, b, transpose = TRUE)
  )
}

#' Weighted sum of the r x r cores: sum_j w_j C_j
#' @keywords internal
#' @noRd
.weighted_core_sum <- function(core, weights) {
  rank_r <- dim(core)[[1L]]
  out <- matrix(0, nrow = rank_r, ncol = rank_r)
  for (j in seq_along(weights)) {
    out <- out + weights[[j]] * core[,, j]
  }
  out
}

#' Rebuild Sigma_j = diag(D_j) + U C_j U^T into a G x G x J array
#' @keywords internal
#' @noRd
.assemble_diag_lowrank_array <- function(diagonal, loadings, core) {
  if (is.null(diagonal) || is.null(loadings) || is.null(core)) {
    stop(
      "`diag_lowrank` requires `diagonal`, `loadings` and `core`.",
      call. = FALSE
    )
  }
  n_genes <- nrow(diagonal)
  n_celltypes <- ncol(diagonal)
  Sigma <- array(0, dim = c(n_genes, n_genes, n_celltypes))
  for (j in seq_len(n_celltypes)) {
    low_rank <- loadings %*% core[,, j] %*% t(loadings)
    Sigma[,, j] <- diag(diagonal[, j], n_genes) + low_rank
  }
  dimnames(Sigma) <- list(
    rownames(loadings),
    rownames(loadings),
    colnames(diagonal)
  )
  Sigma
}

#' Validate a G x G x J covariance array
#' @keywords internal
#' @noRd
.validate_sigma_array <- function(Sigma) {
  if (is.null(Sigma) || length(dim(Sigma)) != 3L) {
    stop("`Sigma` must be a G x G x J array.", call. = FALSE)
  }
  if (dim(Sigma)[[1L]] != dim(Sigma)[[2L]]) {
    stop("`Sigma` slices must be square (G x G).", call. = FALSE)
  }
  invisible(Sigma)
}

#' Normalise and validate block membership
#' @keywords internal
#' @noRd
.validate_covariance_blocks <- function(Sigma, blocks, tol) {
  n_genes <- dim(Sigma)[[1L]]
  if (is.null(blocks)) {
    stop(
      "`structure = \"block\"` requires `blocks`.",
      call. = FALSE
    )
  }
  block_list <- if (is.list(blocks)) {
    lapply(blocks, as.integer)
  } else {
    membership <- as.integer(as.factor(blocks))
    if (length(membership) != n_genes) {
      stop("`blocks` must have length G.", call. = FALSE)
    }
    unname(split(seq_len(n_genes), membership))
  }

  covered <- sort(unlist(block_list))
  if (!identical(covered, seq_len(n_genes))) {
    stop(
      "`blocks` must partition all G genes exactly once.",
      call. = FALSE
    )
  }

  scale <- max(abs(Sigma))
  for (j in seq_len(dim(Sigma)[[3L]])) {
    slice <- Sigma[,, j]
    for (a in seq_along(block_list)) {
      for (bb in seq_along(block_list)) {
        if (a == bb) {
          next
        }
        off <- slice[block_list[[a]], block_list[[bb]], drop = FALSE]
        if (max(abs(off)) > tol * max(scale, 1)) {
          stop(
            "`structure = \"block\"` requires zero between-block ",
            "covariance, but a nonzero cross-block entry was found in ",
            "slice ",
            j,
            ".",
            call. = FALSE
          )
        }
      }
    }
  }
  block_list
}

#' Validate the declared bandwidth against the covariance array
#' @keywords internal
#' @noRd
.validate_covariance_bandwidth <- function(Sigma, bandwidth, tol) {
  if (is.null(bandwidth)) {
    stop(
      "`structure = \"band\"` requires `bandwidth`.",
      call. = FALSE
    )
  }
  bandwidth <- as.integer(bandwidth)
  if (bandwidth < 0L) {
    stop("`bandwidth` must be non-negative.", call. = FALSE)
  }
  n_genes <- dim(Sigma)[[1L]]
  idx <- seq_len(n_genes)
  outside <- abs(outer(idx, idx, "-")) > bandwidth
  scale <- max(abs(Sigma))
  for (j in seq_len(dim(Sigma)[[3L]])) {
    slice <- Sigma[,, j]
    if (max(abs(slice[outside])) > tol * max(scale, 1)) {
      stop(
        "`structure = \"band\"` requires zero entries beyond the ",
        "declared bandwidth, but slice ",
        j,
        " violates it.",
        call. = FALSE
      )
    }
  }
  bandwidth
}

#' Validate diagonal-plus-low-rank factors reconstruct Sigma
#' @keywords internal
#' @noRd
.validate_diag_lowrank <- function(Sigma, diagonal, loadings, core, tol) {
  if (is.null(diagonal) || is.null(loadings) || is.null(core)) {
    stop(
      "`structure = \"diag_lowrank\"` requires `diagonal`, `loadings` ",
      "and `core`.",
      call. = FALSE
    )
  }
  diagonal <- as.matrix(diagonal)
  loadings <- as.matrix(loadings)
  if (nrow(diagonal) != dim(Sigma)[[1L]]) {
    stop("`diagonal` must have G rows.", call. = FALSE)
  }
  if (ncol(diagonal) != dim(Sigma)[[3L]]) {
    stop("`diagonal` must have J columns.", call. = FALSE)
  }
  if (nrow(loadings) != dim(Sigma)[[1L]]) {
    stop("`loadings` must have G rows.", call. = FALSE)
  }
  if (any(diagonal <= 0)) {
    stop("`diagonal` entries must be strictly positive.", call. = FALSE)
  }
  rebuilt <- .assemble_diag_lowrank_array(diagonal, loadings, core)
  scale <- max(abs(Sigma))
  if (max(abs(rebuilt - Sigma)) > tol * max(scale, 1)) {
    stop(
      "`diagonal`, `loadings` and `core` do not reconstruct `Sigma`.",
      call. = FALSE
    )
  }
  list(diagonal = diagonal, loadings = loadings, core = core)
}
