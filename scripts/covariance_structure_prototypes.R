# ---- covariance_structure_prototypes.R -------------------------------------
##
## Prototype the four structure-aware covariance backends that accelerate the
## DeCovarT likelihood when Sigma(p) = sum_j p_j^2 Sigma_j inherits
## exploitable structure. This script is a *prototype*: it verifies each
## backend against the dense Cholesky to machine precision and sketches a
## micro-benchmark, but it does NOT run full deconvolution simulations.
##
## Decision tree (see vignette "generative-model-derivatives",
## section "Structure-aware covariance backends"):
##   dense Cholesky            -> universal default, O(G^3)
##   block Cholesky            -> stochastic block model / disconnected
##                                gene modules (exactly block-diagonal Sigma)
##   banded / sparse Cholesky  -> band / AR, Erdos-Renyi, scale-free,
##                                small-world (sparse covariance; symbolic
##                                ordering computed once, reused per p)
##   Woodbury / Schur          -> shared regulatory programs (pathways) or a
##                                few bridge / housekeeping hub genes
##                                (diagonal-plus-low-rank Sigma)
##
## Run from the package root:
##   Rscript scripts/covariance_structure_prototypes.R

devtools::load_all(".", quiet = TRUE)
stopifnot(requireNamespace("withr", quietly = TRUE))

# Ground-truth dense solve / log-determinant for verification.
dense_reference <- function(p, Sigma) {
  sigma_p <- .compute_global_variance(p, Sigma)
  list(
    logdet = as.numeric(determinant(sigma_p, logarithm = TRUE)$modulus),
    inverse = solve(sigma_p)
  )
}

report_backend <- function(label, cov, p, rhs) {
  ref <- dense_reference(p, cov$Sigma)
  logdet_gap <- abs(sigma_logdet(cov, p) - ref$logdet)
  solve_gap <- max(abs(
    sigma_solve(cov, p, rhs) - as.numeric(ref$inverse %*% rhs)
  ))
  message(
    sprintf(
      "%-14s | logdet gap %.2e | solve gap %.2e",
      label,
      logdet_gap,
      solve_gap
    )
  )
  invisible(list(logdet_gap = logdet_gap, solve_gap = solve_gap))
}

# ===========================================================================
# 1) Block Cholesky -- stochastic block model / disconnected gene modules
# ===========================================================================
# Two gene modules with no between-module covariance in any cell type. The
# assembled Sigma(p) is exactly block-diagonal, so one Cholesky per block
# replaces one dense Cholesky: cost drops by roughly the squared number of
# equal blocks (two equal blocks -> ~4x fewer factorisation flops).

make_spd <- function(n, seed) {
  withr::with_seed(seed, {
    a <- matrix(stats::rnorm(n * n), n, n)
    crossprod(a) + diag(n)
  })
}

block_sizes <- c(3L, 3L)
n_genes_block <- sum(block_sizes)
membership <- rep(seq_along(block_sizes), block_sizes)
Sigma_block <- array(0, dim = c(n_genes_block, n_genes_block, 2L))
for (j in 1:2) {
  offset <- 0L
  for (b in seq_along(block_sizes)) {
    idx <- offset + seq_len(block_sizes[b])
    Sigma_block[idx, idx, j] <- make_spd(block_sizes[b], 100L * j + b)
    offset <- offset + block_sizes[b]
  }
}

cov_block <- new_decovart_covariance(
  Sigma_block,
  structure = "block",
  blocks = membership
)
report_backend(
  "block",
  cov_block,
  p = c(0.6, 0.4),
  rhs = seq_len(n_genes_block)
)

# ===========================================================================
# 2) Band matrix -- AR(1)-like ordered local dependence
# ===========================================================================
# Only neighbouring genes (|i - j| <= b) covary. A banded Cholesky costs
# O(G b^2) instead of O(G^3). Routed through the sparse backend.

build_band <- function(n, bandwidth, seed) {
  m <- diag(n) * 3
  withr::with_seed(seed, {
    for (i in seq_len(n - 1L)) {
      for (k in seq_len(min(bandwidth, n - i))) {
        val <- stats::runif(1, -0.3, 0.3)
        m[i, i + k] <- val
        m[i + k, i] <- val
      }
    }
  })
  m
}
n_genes_band <- 8L
Sigma_band <- array(0, dim = c(n_genes_band, n_genes_band, 2L))
Sigma_band[,, 1] <- build_band(n_genes_band, 1L, 201L)
Sigma_band[,, 2] <- build_band(n_genes_band, 1L, 202L)

cov_band <- new_decovart_covariance(
  Sigma_band,
  structure = "band",
  bandwidth = 1L
)
report_backend(
  "band",
  cov_band,
  p = c(0.5, 0.5),
  rhs = seq_len(n_genes_band)
)

# ===========================================================================
# 3) Sparse matrix -- Erdos-Renyi sparsity + symbolic factorisation
# ===========================================================================
# A sparse covariance with a fixed nonzero pattern across p. The
# fill-reducing ordering (symbolic factorisation) is computed once and reused
# for every numeric refactorisation as p changes during optimisation.
# NOTE: a sparse gene-network PRECISION does not imply a sparse COVARIANCE;
# this backend accelerates the likelihood only when Sigma(p) itself is
# sparse.

build_er_sparse <- function(n, seed) {
  m <- diag(n) * 5
  withr::with_seed(seed, {
    n_edges <- n
    for (e in seq_len(n_edges)) {
      ij <- sample.int(n, 2L)
      val <- stats::runif(1, -0.4, 0.4)
      m[ij[1], ij[2]] <- val
      m[ij[2], ij[1]] <- val
    }
  })
  m
}
n_genes_sparse <- 12L
Sigma_sparse <- array(0, dim = c(n_genes_sparse, n_genes_sparse, 2L))
Sigma_sparse[,, 1] <- build_er_sparse(n_genes_sparse, 301L)
Sigma_sparse[,, 2] <- build_er_sparse(n_genes_sparse, 302L)

cov_sparse <- new_decovart_covariance(Sigma_sparse, structure = "sparse")
# First solve builds the ordering; the second reuses it (same object).
report_backend(
  "sparse (p1)",
  cov_sparse,
  p = c(0.6, 0.4),
  rhs = seq_len(n_genes_sparse)
)
report_backend(
  "sparse (p2)",
  cov_sparse,
  p = c(0.35, 0.65),
  rhs = seq_len(n_genes_sparse)
)

# ===========================================================================
# 4) Intrinsic low rank -- shared regulatory programs (Schur / Woodbury)
# ===========================================================================
# Genes interact through a small number r << G of shared latent programs
# (pathways), plus gene-specific noise: Sigma_j = D_j + U C_j U^T. This is the
# factor-analysis form Sigma = W W^T + Psi. Woodbury moves the O(G^3) solve to
# the low rank r; the matrix-determinant lemma does the same for the
# log-determinant. The "few bridge / housekeeping hub genes" reading is the
# arrow-matrix special case handled by the same Schur-complement algebra.

n_genes_lr <- 20L
rank_r <- 3L
n_ct <- 2L
diagonal_lr <- matrix(
  withr::with_seed(401L, stats::runif(n_genes_lr * n_ct, 0.5, 1.5)),
  nrow = n_genes_lr
)
loadings_lr <- matrix(
  withr::with_seed(402L, stats::rnorm(n_genes_lr * rank_r)),
  nrow = n_genes_lr
)
core_lr <- array(0, dim = c(rank_r, rank_r, n_ct))
core_lr[,, 1] <- make_spd(rank_r, 403L)
core_lr[,, 2] <- make_spd(rank_r, 404L)

cov_lr <- new_decovart_covariance(
  structure = "diag_lowrank",
  diagonal = diagonal_lr,
  loadings = loadings_lr,
  core = core_lr
)
report_backend(
  "diag_lowrank",
  cov_lr,
  p = c(0.55, 0.45),
  rhs = seq_len(n_genes_lr)
)

# ---------------------------------------------------------------------------
# Optional micro-benchmark scaffold (prototype only; not a full simulation).
# Uncomment to compare a single factorise+solve for dense vs a structured
# backend at a larger G. Keep G modest so the script stays fast.
# ---------------------------------------------------------------------------
# if (requireNamespace("bench", quietly = TRUE)) {
#   p_bench <- c(0.5, 0.5)
#   rhs_bench <- seq_len(n_genes_lr)
#   print(bench::mark(
#     dense = {
#       s <- .compute_global_variance(p_bench, cov_lr$Sigma)
#       solve(s, rhs_bench)
#     },
#     woodbury = sigma_solve(cov_lr, p_bench, rhs_bench),
#     check = FALSE
#   ))
# }

message("All covariance-structure prototypes verified against dense.")
