#' Sobol points in the unit cube, with a uniform fallback
#'
#' @keywords internal
#' @noRd
.sobol_unit_cube <- function(n, d, seed = NULL) {
  n <- as.integer(n)
  d <- as.integer(d)
  if (n < 1L || d < 1L) {
    stop("`n` and `d` must be positive integers.", call. = FALSE)
  }
  if (requireNamespace("qrng", quietly = TRUE)) {
    args <- list(n = n, d = d, randomize = "digital.shift")
    if (!is.null(seed)) {
      args$seed <- as.integer(seed)
    }
    u <- do.call(qrng::sobol, args)
    if (is.null(dim(u))) {
      u <- matrix(u, ncol = 1L)
    }
    return(u)
  }
  .ui_warn(
    "{.pkg qrng} is not installed; Sobol points fall back to",
    " {.fn stats::runif}."
  )
  if (!is.null(seed)) {
    set.seed(seed)
  }
  matrix(stats::runif(n * d), nrow = n, ncol = d)
}

#' Inverse-transform Gaussian draws from unit-cube points
#'
#' @keywords internal
#' @noRd
.mvnorm_from_unit <- function(u, mu, chol_upper) {
  eps <- 1e-12
  u <- pmin(pmax(u, eps), 1 - eps)
  z <- stats::qnorm(u)
  mu <- as.numeric(mu)
  sweep(z %*% chol_upper, 2L, mu, "+")
}

#' Log-density of rows under N(mu, R^T R)
#'
#' @keywords internal
#' @noRd
.log_dmvnorm_chol <- function(x, mu, chol_upper, log_det) {
  mu <- as.numeric(mu)
  g <- length(mu)
  z <- backsolve(
    chol_upper,
    t(x) - mu,
    transpose = TRUE
  )
  -0.5 * (g * log(2 * pi) + log_det + colSums(z * z))
}

#' Cholesky factor and log-determinant of one covariance
#'
#' @keywords internal
#' @noRd
.chol_logdet <- function(sigma) {
  sigma <- as.matrix(sigma)
  chol_upper <- chol(sigma)
  list(
    chol = chol_upper,
    log_det = 2 * sum(log(diag(chol_upper)))
  )
}

#' MixSim-style pairwise overlap via stratified Sobol Monte Carlo
#'
#' Estimates the MixSim Omega map
#' \eqn{\Omega_{j\ell}=\Pr_{X\sim f_j}(\pi_\ell f_\ell(X)>\pi_j f_j(X))}
#' by drawing `n_mc` quasi-Monte Carlo points **from each component**
#' (inverse-transform sampling of a Sobol sequence through the
#' Cholesky factor) and comparing **log-densities**. The average
#' pairwise overlap `BarOmega` is the unweighted mean of
#' \eqn{\Omega_{j\ell}+\Omega_{\ell j}} over \eqn{j<\ell}, matching
#' [MixSim::overlap()] / [FSDA MixSim](https://rosa.unipr.it/FSDA/MixSim.html).
#'
#' For two equal-weight components this pairwise overlap equals the
#' histogram similarity \eqn{\int\min(f_j,f_\ell)} and therefore
#' \eqn{1-\mathrm{TV}(f_j,f_\ell)}
#' \insertCite{nielsenGuaranteedDeterministicBounds2018}{DeCovarT}.
#'
#' @param true_theta List with `p`, `mu`, `sigma` as in
#'   [check_true_theta()].
#' @param n_mc Draws **per component** (default `10000`).
#' @param seed Optional RNG seed (passed to Sobol randomisation and to
#'   the uniform fallback).
#' @param J Optional number of components.
#'
#' @return A list with `BarOmega`, `MaxOmega`, and `OmegaMap`.
#' @export
#' @references
#' \insertAllCited{}
#' @seealso [compute_average_overlap()], [compute_average_riemannian()]
#' @examples
#' set.seed(1)
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' overlap_gaussian_mc(theta, n_mc = 2000L)$BarOmega
overlap_gaussian_mc <- function(
  true_theta,
  n_mc = 10000L,
  seed = NULL,
  J = NULL
) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = TRUE,
    J = J,
    second_moment = "sigma"
  )
  n_mc <- as.integer(n_mc)
  if (length(n_mc) != 1L || is.na(n_mc) || n_mc < 1L) {
    stop("`n_mc` must be a positive integer.", call. = FALSE)
  }
  p <- as.numeric(theta$p)
  mu <- theta$mu
  sigma <- theta$sigma
  n_celltypes <- theta$J
  n_genes <- theta$G
  log_pi <- log(pmax(p, .Machine$double.xmin))

  factors <- vector("list", n_celltypes)
  for (j in seq_len(n_celltypes)) {
    factors[[j]] <- .chol_logdet(sigma[,, j])
  }

  draws <- vector("list", n_celltypes)
  for (j in seq_len(n_celltypes)) {
    stream_seed <- if (is.null(seed)) {
      NULL
    } else {
      as.integer(seed) + 1000L * j
    }
    u <- .sobol_unit_cube(n_mc, n_genes, seed = stream_seed)
    draws[[j]] <- .mvnorm_from_unit(
      u,
      mu[, j],
      factors[[j]]$chol
    )
  }

  omega_map <- diag(n_celltypes)
  n_pairs <- as.integer(n_celltypes * (n_celltypes - 1L) / 2L)
  pair_omega <- numeric(n_pairs)
  pair_idx <- 0L
  for (j in seq_len(n_celltypes - 1L)) {
    for (ell in seq.int(j + 1L, n_celltypes)) {
      lj_on_j <- .log_dmvnorm_chol(
        draws[[j]],
        mu[, j],
        factors[[j]]$chol,
        factors[[j]]$log_det
      )
      le_on_j <- .log_dmvnorm_chol(
        draws[[j]],
        mu[, ell],
        factors[[ell]]$chol,
        factors[[ell]]$log_det
      )
      le_on_e <- .log_dmvnorm_chol(
        draws[[ell]],
        mu[, ell],
        factors[[ell]]$chol,
        factors[[ell]]$log_det
      )
      lj_on_e <- .log_dmvnorm_chol(
        draws[[ell]],
        mu[, j],
        factors[[j]]$chol,
        factors[[j]]$log_det
      )
      w_ell_j <- mean(log_pi[[ell]] + le_on_j > log_pi[[j]] + lj_on_j)
      w_j_ell <- mean(log_pi[[j]] + lj_on_e > log_pi[[ell]] + le_on_e)
      omega_map[j, ell] <- w_ell_j
      omega_map[ell, j] <- w_j_ell
      pair_idx <- pair_idx + 1L
      pair_omega[[pair_idx]] <- w_ell_j + w_j_ell
    }
  }
  bar <- if (n_pairs < 1L) {
    NA_real_
  } else {
    mean(pair_omega)
  }
  list(
    BarOmega = bar,
    MaxOmega = if (n_pairs < 1L) {
      NA_real_
    } else {
      max(pair_omega)
    },
    OmegaMap = omega_map
  )
}

#' Affine-invariant Riemannian (AIRM) distance between two SPD matrices
#'
#' AIRM is the **affine-invariant Riemannian metric** on the cone of
#' symmetric positive-definite matrices:
#' \eqn{d_R(A,B)=\lVert\log(A^{-1/2}BA^{-1/2})\rVert_F}.
#' Do **not** use the Frobenius
#' \eqn{\lVert A-B\rVert_F}: that Euclidean chord leaves the cone,
#' is not inversion-invariant, and suffers the swelling effect.
#'
#' @param a,b Symmetric positive-definite matrices of equal size.
#' @return Non-negative scalar.
#' @export
#' @seealso [compute_average_riemannian()]
spd_affine_invariant_distance <- function(a, b) {
  a <- as.matrix(a)
  b <- as.matrix(b)
  if (!is.numeric(a) || !is.numeric(b) || anyNA(a) || anyNA(b)) {
    stop("`a` and `b` must be numeric matrices without NA.", call. = FALSE)
  }
  if (!identical(dim(a), dim(b)) || nrow(a) != ncol(a)) {
    stop("`a` and `b` must be square matrices of equal size.", call. = FALSE)
  }
  r <- chol(a)
  # C = R^{-T} B R^{-1} with A = R^T R (upper Cholesky).
  tmp <- backsolve(r, b, transpose = TRUE)
  c_mat <- backsolve(r, t(tmp), transpose = TRUE)
  ev <- eigen(c_mat, symmetric = TRUE, only.values = TRUE)$values
  ev <- pmax(ev, .Machine$double.xmin)
  sqrt(sum(log(ev)^2))
}

#' Mean pairwise AIRM distance of component covariances
#'
#' @inheritParams compute_average_jeffreys
#' @return Scalar average of \eqn{d_R(\Sigma_j,\Sigma_\ell)} over
#'   \eqn{j<\ell}.
#' @export
#' @seealso [spd_affine_invariant_distance()]
compute_average_riemannian <- function(true_theta, J = NULL) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = FALSE,
    J = J,
    second_moment = "sigma"
  )
  sigma <- theta$sigma
  n_celltypes <- theta$J
  acc <- 0
  n_pairs <- 0L
  for (j in seq_len(n_celltypes - 1L)) {
    for (ell in seq.int(j + 1L, n_celltypes)) {
      acc <- acc +
        spd_affine_invariant_distance(sigma[,, j], sigma[,, ell])
      n_pairs <- n_pairs + 1L
    }
  }
  if (n_pairs < 1L) {
    return(NA_real_)
  }
  acc / n_pairs
}

#' Scale component covariances, keeping precision zeros
#'
#' Replaces each \eqn{\Sigma_j} by \eqn{s\Sigma_j}. Then
#' \eqn{\Omega_j(s)=\Omega_j/s}, so the **support** of the precision
#' (structural zeros from the graph) is unchanged.
#'
#' @param sigma \eqn{G\times G\times J} covariance array.
#' @param scale Positive scalar \eqn{s}.
#' @return Scaled array.
#' @export
#' @seealso [scale_covariances_to_overlap()]
scale_covariance_array <- function(sigma, scale) {
  scale <- as.numeric(scale)
  if (length(scale) != 1L || !is.finite(scale) || scale <= 0) {
    stop("`scale` must be a single positive number.", call. = FALSE)
  }
  .assert_ggj_array(sigma, "sigma")
  sigma * scale
}

#' Invert a covariance array to precisions
#'
#' @keywords internal
#' @noRd
.precision_array_from_sigma <- function(sigma) {
  d <- dim(sigma)
  theta <- array(NA_real_, dim = d)
  for (j in seq_len(d[[3L]])) {
    theta[,, j] <- solve(sigma[,, j])
  }
  theta
}

#' Calibrate a global covariance scale to a target MixSim BarOmega
#'
#' Binary-searches \eqn{s>0} so that
#' [compute_average_overlap()] of \eqn{\{s\Sigma_j\}} matches
#' `target`. Larger \eqn{s} inflates the Gaussians and **increases**
#' average overlap while leaving graph zeros intact. This is the
#' graph-constrained analogue of MixSim / FSDA's average-overlap
#' generator, which cannot take a precision support as input.
#'
#' @param mu \eqn{G\times J} mean signature.
#' @param sigma Base \eqn{G\times G\times J} covariances.
#' @param p Simplex weights (length \eqn{J}).
#' @param target Target average pairwise overlap in \eqn{(0,1)}.
#' @param n_mc Monte Carlo size forwarded to
#'   [overlap_gaussian_mc()] when \eqn{G\ge 4}.
#' @param s_range Length-2 search interval for \eqn{s}.
#' @param tol Absolute overlap tolerance.
#' @param seed Optional seed for overlap evaluations.
#' @param verbose If `TRUE`, print the MixSim-to-MC switch.
#'
#' @return A list with `sigma`, `Theta`, `scale`, `baromega`, and
#'   `target`.
#' @export
scale_covariances_to_overlap <- function(
  mu,
  sigma,
  p,
  target,
  n_mc = 4000L,
  s_range = c(1e-3, 1e3),
  tol = 0.005,
  seed = NULL,
  verbose = FALSE
) {
  target <- as.numeric(target)
  if (
    length(target) != 1L || !is.finite(target) || target <= 0 || target >= 1
  ) {
    stop("`target` must lie in (0, 1).", call. = FALSE)
  }
  s_range <- as.numeric(s_range)
  if (
    length(s_range) != 2L ||
      any(s_range <= 0) ||
      s_range[[1L]] >= s_range[[2L]]
  ) {
    stop(
      "`s_range` must be two increasing positive numbers.",
      call. = FALSE
    )
  }
  overlap_at <- function(s) {
    th <- list(
      p = p,
      mu = mu,
      sigma = scale_covariance_array(sigma, s)
    )
    compute_average_overlap(
      th,
      n_mc = n_mc,
      seed = seed,
      verbose = verbose
    )
  }
  ov_lo <- overlap_at(s_range[[1L]])
  ov_hi <- overlap_at(s_range[[2L]])
  if (!is.finite(ov_lo) || !is.finite(ov_hi)) {
    stop("Overlap evaluations at `s_range` were not finite.", call. = FALSE)
  }
  if (target <= min(ov_lo, ov_hi) || target >= max(ov_lo, ov_hi)) {
    s_star <- if (abs(ov_lo - target) <= abs(ov_hi - target)) {
      s_range[[1L]]
    } else {
      s_range[[2L]]
    }
    .ui_warn(c(
      "Target overlap {.val {target}} is outside the range",
      "[{.val {round(ov_lo, 4)}}, {.val {round(ov_hi, 4)}}]",
      "on `s_range`; using the nearer endpoint."
    ))
  } else {
    s_star <- stats::uniroot(
      function(s) overlap_at(s) - target,
      interval = s_range,
      tol = max(tol / 4, 1e-4),
      maxiter = 40L
    )$root
  }
  sigma_s <- scale_covariance_array(sigma, s_star)
  bar <- overlap_at(s_star)
  list(
    sigma = sigma_s,
    Theta = .precision_array_from_sigma(sigma_s),
    scale = s_star,
    baromega = bar,
    target = target
  )
}
