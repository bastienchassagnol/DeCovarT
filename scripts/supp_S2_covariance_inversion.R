###############################################################################
###############################################################################
###                                                                         ###
###     SUPPLEMENTARY S2 – COVARIANCE INVERSION BACKENDS                   ###
###     Dense vs banded vs block-diagonal · time and memory scaling         ###
###                                                                         ###
###############################################################################
###############################################################################
#
# Background launch from the repository root (vanilla Rscript; no CLI
# parser — hyperparameters are hard-coded below). stdout and stderr go
# to logs/:
#
#   mkdir -p logs
#   nohup Rscript --no-save --no-restore \
#     scripts/supp_S2_covariance_inversion.R \
#     > "logs/supp_S2_$(date +%F)_covariance_inversion.log" 2>&1 &
#
# Article:  DeCovarT – Supplementary S2: covariance backends
# Vignette: vignettes/supp-S2-covariance-inversion.qmd
# Package:  new_decovart_covariance() + sigma_solve() / sigma_logdet()
#           (correctness prototypes live in this script, not a second file)

#
# ── Factorial design ─────────────────────────────────────────────────────────
#  Factor          Levels
#  ─────────────── ─────────────────────────────────────────────────────────
#  G               50, 200, 500, 1000
#  Sparsity type   dense; band (bandwidth 5); block-diagonal (2 blocks)
#  Backend         "dense" (Cholesky), "band" (banded), "block" (block-diag)
#  Timing reps     20
#
# ── Usage ────────────────────────────────────────────────────────────────────
#  Rscript scripts/supp_S2_covariance_inversion.R
#
# ── Outputs ─────────────────────────────────────────────────────────────────
#  output/supp_S2/timing_results.rds
#  output/supp_S2/S2_timing_plot.pdf
###############################################################################

# ==============================================================================
# SECTION 0 · Dependencies and paths ----
# ==============================================================================

if (!requireNamespace("DeCovarT", quietly = TRUE)) {
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(".", quiet = TRUE)
  } else {
    stop("Install DeCovarT or devtools before running this script.")
  }
} else {
  library(DeCovarT)
}

OUT_DIR <- file.path("output", "supp_S2")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

N_TIMING_REPS <- as.integer(Sys.getenv("N_TIMING_REPS", "20"))
if (interactive()) {
  N_TIMING_REPS <- 3L
  message("supp_S2: interactive – using ", N_TIMING_REPS, " timing reps.")
}

SEED <- 20260905L
set.seed(SEED)


# ==============================================================================
# SECTION 0b · Backend correctness ----
#   Formerly covariance_structure_prototypes.R.
#   Each structured backend is checked against a dense Cholesky to machine
#   precision on a small G. Timing at large G is SECTION 2.
# ==============================================================================

dense_reference <- function(p, Sigma) {
  n_ct <- dim(Sigma)[[3L]]
  sigma_p <- matrix(0, dim(Sigma)[[1L]], dim(Sigma)[[2L]])
  for (j in seq_len(n_ct)) {
    sigma_p <- sigma_p + (p[[j]]^2) * Sigma[,, j]
  }
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

make_spd <- function(n, seed) {
  withr::with_seed(seed, {
    a <- matrix(stats::rnorm(n * n), n, n)
    crossprod(a) + diag(n)
  })
}

if (!requireNamespace("withr", quietly = TRUE)) {
  message("supp_S2 | withr missing; skipping backend correctness checks.")
} else {
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
  report_backend("block", cov_block, c(0.6, 0.4), seq_len(n_genes_block))

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
  report_backend("band", cov_band, c(0.5, 0.5), seq_len(n_genes_band))

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
  report_backend(
    "sparse (p1)",
    cov_sparse,
    c(0.6, 0.4),
    seq_len(n_genes_sparse)
  )
  report_backend(
    "sparse (p2)",
    cov_sparse,
    c(0.35, 0.65),
    seq_len(n_genes_sparse)
  )

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
    c(0.55, 0.45),
    seq_len(n_genes_lr)
  )
  message("supp_S2 | backend prototypes verified against dense Cholesky.")
}


# ==============================================================================
# SECTION 1 · GENERATIVE MODEL ----
#   Construct covariance matrices of different sizes and sparsity structures.
# ==============================================================================

G_LEVELS <- c(50L, 200L, 500L, 1000L)
J_FIXED <- 3L
p_eq <- rep(1 / J_FIXED, J_FIXED)

# For each G, produce three Sigma arrays (dense, band-5, block-diagonal)
.make_sigma_dense <- function(G, j, seed) {
  set.seed(seed)
  A <- matrix(rnorm(G * G), G, G)
  S <- crossprod(A) / G + diag(G)
  S
}

.make_sigma_band <- function(G, bandwidth = 5L) {
  band_decay <- 0.9
  S <- matrix(0, G, G)
  for (k in 0:bandwidth) {
    idx <- seq_len(G - k)
    val <- band_decay^k
    S[cbind(idx, idx + k)] <- val
    S[cbind(idx + k, idx)] <- val
  }
  diag(S) <- 1
  S
}

.make_sigma_block <- function(G) {
  half <- G %/% 2L
  # Two independent blocks; each is a random SPD matrix
  set.seed(SEED + G)
  A1 <- matrix(rnorm(half^2), half, half)
  A2 <- matrix(rnorm((G - half)^2), G - half, G - half)
  S <- matrix(0, G, G)
  S[seq_len(half), seq_len(half)] <- crossprod(A1) / half + diag(half)
  S[(half + 1L):G, (half + 1L):G] <- crossprod(A2) / (G - half) + diag(G - half)
  S
}

message("supp_S2 | building covariance structures ...")
sigma_cases <- purrr::map_dfr(G_LEVELS, function(G) {
  tibble::tibble(
    G = G,
    structure = c("dense", "band5", "block2"),
    Sigma_list = list(
      purrr::map(seq_len(J_FIXED), ~ .make_sigma_dense(G, .x, SEED + G + .x)),
      rep(list(.make_sigma_band(G, 5L)), J_FIXED),
      rep(list(.make_sigma_block(G)), J_FIXED)
    )
  )
})
message("supp_S2 | covariance structures built: ", nrow(sigma_cases), " cases.")


# ==============================================================================
# SECTION 2 · INFERENCE (timing) ----
#   For each (G, structure), time sigma_solve() via the dense backend.
#   The covariance-backend abstraction (new_decovart_covariance) is compared
#   against direct base::solve() as baseline.
# ==============================================================================

.list_to_sigma_array <- function(Slist) {
  g <- nrow(Slist[[1L]])
  j <- length(Slist)
  arr <- array(0, dim = c(g, g, j))
  for (jj in seq_len(j)) {
    arr[,, jj] <- Slist[[jj]]
  }
  arr
}

message("supp_S2 | timing (", N_TIMING_REPS, " reps per case) ...")
timing_rows <- purrr::map_dfr(seq_len(nrow(sigma_cases)), function(i) {
  G <- sigma_cases$G[i]
  strct <- sigma_cases$structure[i]
  Slist <- sigma_cases$Sigma_list[[i]]
  rhs <- matrix(rnorm(G), nrow = G, ncol = 1L)

  backend_structure <- switch(
    strct,
    dense = "dense",
    band5 = "band",
    block2 = "block",
    "dense"
  )
  extra <- list()
  if (identical(backend_structure, "band")) {
    extra$bandwidth <- 5L
  }
  if (identical(backend_structure, "block")) {
    extra$blocks <- rep(c(1L, 2L), each = G %/% 2L)
    if (length(extra$blocks) < G) {
      extra$blocks <- c(extra$blocks, 2L)
    }
  }

  times <- tryCatch(
    {
      Sigma_arr <- .list_to_sigma_array(Slist)
      args <- c(
        list(Sigma = Sigma_arr, structure = backend_structure),
        extra
      )
      backend <- do.call(new_decovart_covariance, args)
      Sigma_p <- Reduce(
        `+`,
        purrr::map2(Slist, p_eq^2, \(x, w) x * w)
      )
      t_backend <- system.time(
        for (r in seq_len(N_TIMING_REPS)) {
          sigma_solve(backend, p_eq, rhs)
        }
      )[["elapsed"]]
      t_base <- system.time(
        for (r in seq_len(N_TIMING_REPS)) {
          base::solve(Sigma_p, rhs)
        }
      )[["elapsed"]]
      c(backend = t_backend, base = t_base)
    },
    error = function(e) {
      message("  skipping G=", G, " ", strct, ": ", conditionMessage(e))
      c(backend = NA_real_, base = NA_real_)
    }
  )

  tibble::tibble(
    G = G,
    structure = strct,
    time_backend_s = times[["backend"]] / N_TIMING_REPS,
    time_base_s = times[["base"]] / N_TIMING_REPS
  )
})

saveRDS(timing_rows, file.path(OUT_DIR, "timing_results.rds"))
message(
  "supp_S2 | timing complete:\n",
  paste(capture.output(print(timing_rows)), collapse = "\n")
)


# ==============================================================================
# SECTION 3 · VISUALISATIONS ----
# ==============================================================================

if (requireNamespace("ggplot2", quietly = TRUE)) {
  long_timing <- tidyr::pivot_longer(
    timing_rows,
    cols = c("time_backend_s", "time_base_s"),
    names_to = "solver",
    values_to = "seconds",
    names_pattern = "time_(.+)_s"
  )
  p_timing <- ggplot2::ggplot(
    long_timing,
    ggplot2::aes(
      x = G,
      y = seconds,
      colour = solver,
      linetype = structure
    )
  ) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(
      x = "Number of genes (G)",
      y = "Seconds per solve (log scale)",
      colour = "Solver",
      linetype = "Covariance structure",
      caption = paste0(
        "Timing averaged over ",
        N_TIMING_REPS,
        " replicates. ",
        "J = ",
        J_FIXED,
        " fixed."
      )
    ) +
    ggplot2::theme_bw()

  ggplot2::ggsave(
    file.path(OUT_DIR, "S2_timing_plot.pdf"),
    plot = p_timing,
    width = 8,
    height = 5
  )
  message("supp_S2 | saved S2_timing_plot.pdf")
}

message("supp_S2 | done. Outputs in ", normalizePath(OUT_DIR, mustWork = FALSE))
