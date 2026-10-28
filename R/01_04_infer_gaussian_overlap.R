#' Parse and standardise \eqn{\theta} (internal)
#'
#' `true_theta` is a named list with `p` (length \eqn{J}, or a
#' \eqn{J\times N} / \eqn{N\times J} matrix of sample-wise proportions, averaged
#' over samples), `mu` (\eqn{G\times J}, or \eqn{J\times G}), and `sigma`
#' and/or `Theta` (\eqn{G\times G\times J} covariances / precisions).
#'
#' @return Named list with `p` (length-\eqn{J} vector or `NULL`),
#'   `mu` (\eqn{G\times J}), `sigma` and/or `Theta` (\eqn{G\times G\times J}),
#'   `G`, and `J`.
#' @keywords internal
#' @noRd
.parse_true_theta <- function(
  true_theta,
  require_p = TRUE,
  J = NULL,
  second_moment = c("either", "sigma", "Theta")
) {
  second_moment <- .match_arg_case_insensitive(
    second_moment,
    c("either", "sigma", "Theta")
  )
  if (!is.list(true_theta)) {
    stop("`true_theta` must be a list.", call. = FALSE)
  }
  if (is.null(true_theta$mu)) {
    stop("`true_theta` must contain `mu`.", call. = FALSE)
  }

  # --- Second moments: require Sigma and/or precision Theta -------------------
  # Each must be a G x G x J array (one PD matrix per cell type). When both are
  # supplied they must share the same (G, J).
  has_sigma <- !is.null(true_theta$sigma)
  has_theta <- !is.null(true_theta$Theta)
  if (identical(second_moment, "sigma") && !has_sigma) {
    stop("`true_theta` must contain `sigma`.", call. = FALSE)
  }
  if (identical(second_moment, "Theta") && !has_theta) {
    stop("`true_theta` must contain `Theta`.", call. = FALSE)
  }
  if (identical(second_moment, "either") && !has_sigma && !has_theta) {
    stop(
      "`true_theta` must contain `sigma` and/or `Theta`.",
      call. = FALSE
    )
  }

  # Helper: enforce cubic layout G x G x J and return (G, J).
  # Shared implementation lives in utils.R as `.assert_ggj_array()`.
  .check_ggj_array <- function(arr, nm) {
    .assert_ggj_array(arr, paste0("true_theta$", nm))
  }

  dims_list <- list()
  if (has_sigma) {
    dims_list$sigma <- .check_ggj_array(true_theta$sigma, "sigma")
  }
  if (has_theta) {
    dims_list$Theta <- .check_ggj_array(true_theta$Theta, "Theta")
  }
  if (length(dims_list) == 2L) {
    if (
      !identical(dims_list$sigma$G, dims_list$Theta$G) ||
        !identical(dims_list$sigma$J, dims_list$Theta$J)
    ) {
      stop("`sigma` and `Theta` must share the same G x G x J dims.")
    }
  }
  G <- dims_list[[1L]]$G
  n_celltypes <- dims_list[[1L]]$J

  # --- Mean signature mu: genes x cell types ---------------------------------
  # Preferred storage is G x J; also accept MixSim-style J x G and transpose.
  mu <- as.matrix(true_theta$mu)
  if (!is.numeric(mu) || anyNA(mu)) {
    stop("`true_theta$mu` must be numeric without missing values.")
  }
  if (nrow(mu) == G && ncol(mu) == n_celltypes) {
    mu_gj <- mu
  } else if (nrow(mu) == n_celltypes && ncol(mu) == G) {
    mu_gj <- t(mu)
  } else {
    stop("`true_theta$mu` must be G x J (or J x G).")
  }

  # --- Mixing proportions p: unit simplex ------------------------------------
  # Accept a single mixture weight vector (length J) or sample-wise matrices
  # (J x N or N x J). Matrix forms are reduced to length-J by averaging so that
  # MixSim / Jeffreys summaries still see one weight per cell type. Every
  # supplied simplex (vector, or each column/row of a matrix) must be
  # non-negative and sum to 1.
  p_vec <- NULL
  if (is.null(true_theta$p)) {
    if (isTRUE(require_p)) {
      stop("`true_theta` must contain `p`.", call. = FALSE)
    }
  } else {
    p_raw <- true_theta$p
    if (is.null(dim(p_raw))) {
      # Length-J vector on the simplex (shared mixture weights).
      p_vec <- as.numeric(p_raw)
      if (length(p_vec) != n_celltypes) {
        stop(
          "`true_theta$p` must have length J (= third dim of sigma/Theta).",
          call. = FALSE
        )
      }
    } else {
      p_mat <- as.matrix(p_raw)
      if (nrow(p_mat) == n_celltypes) {
        # J x N: each column is a sample-wise cellular ratio on the simplex;
        # average across samples for mixture-level summaries.
        if (
          any(abs(colSums(p_mat) - 1) > 100 * .Machine$double.eps) ||
            any(p_mat < 0) ||
            anyNA(p_mat)
        ) {
          stop(
            "`true_theta$p` as J x N must have non-negative columns ",
            "summing to 1.",
            call. = FALSE
          )
        }
        p_vec <- rowMeans(p_mat)
      } else if (ncol(p_mat) == n_celltypes) {
        # N x J: each row is a sample-wise ratio; average across samples.
        if (
          any(abs(rowSums(p_mat) - 1) > 100 * .Machine$double.eps) ||
            any(p_mat < 0) ||
            anyNA(p_mat)
        ) {
          stop(
            "`true_theta$p` as N x J must have non-negative rows ",
            "summing to 1.",
            call. = FALSE
          )
        }
        p_vec <- colMeans(p_mat)
      } else {
        stop(
          "`true_theta$p` must be length J, J x N, or N x J.",
          call. = FALSE
        )
      }
    }
    # Final simplex check on the (possibly averaged) length-J weight vector.
    if (
      any(p_vec < 0) ||
        anyNA(p_vec) ||
        abs(sum(p_vec) - 1) > 100 * .Machine$double.eps
    ) {
      stop(
        "`true_theta$p` must be non-negative and sum to 1 ",
        "(after averaging over samples if matrix).",
        call. = FALSE
      )
    }
  }

  # --- Cell-type count J -----------------------------------------------------
  # Infer from the third dim of sigma/Theta, or verify a caller-supplied J.
  if (is.null(J)) {
    J <- n_celltypes
  } else {
    J <- as.integer(J)
    if (length(J) != 1L || is.na(J) || J != n_celltypes) {
      stop(
        "`J` must match the third dimension of sigma/Theta.",
        call. = FALSE
      )
    }
  }
  if (J < 2L) {
    stop("`J` must be at least 2.", call. = FALSE)
  }

  # Standardised fields for downstream metrics (mu always G x J).
  list(
    p = p_vec,
    mu = mu_gj,
    sigma = if (has_sigma) true_theta$sigma else NULL,
    Theta = if (has_theta) true_theta$Theta else NULL,
    G = G,
    J = J
  )
}

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
  cli::cli_inform(
    "{.pkg qrng} is not installed; Sobol points use {.fn stats::runif}.",
    .frequency = "once",
    .frequency_id = "DeCovarT_sobol_qrng"
  )
  .with_seed(seed, matrix(stats::runif(n * d), nrow = n, ncol = d))
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
#' @param true_theta List with `p` (length \eqn{J}), `mu` (\eqn{G\times J})
#'   and `sigma` (\eqn{G\times G\times J}).
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

#' Mean pairwise AIRM distance of component covariances
#'
#' Average over \eqn{j<\ell} of the **affine-invariant Riemannian metric**
#' on the cone of symmetric positive-definite matrices,
#' \eqn{d_R(\Sigma_j,\Sigma_\ell)=\lVert\log(\Sigma_j^{-1/2}\Sigma_\ell
#' \Sigma_j^{-1/2})\rVert_F}, computed from the eigenvalues of
#' \eqn{R^{-\mathsf{T}}\Sigma_\ell R^{-1}} with \eqn{\Sigma_j=R^{\mathsf{T}}R}.
#' The Frobenius chord \eqn{\lVert\Sigma_j-\Sigma_\ell\rVert_F} is not used:
#' it leaves the cone, is not inversion-invariant, and suffers the swelling
#' effect. With \eqn{J=2} the average is the pairwise distance.
#'
#' @inheritParams compute_average_jeffreys
#' @return Scalar average of \eqn{d_R(\Sigma_j,\Sigma_\ell)} over
#'   \eqn{j<\ell}.
#' @export
#' @examples
#' a <- diag(2)
#' b <- matrix(c(2, 0.5, 0.5, 1), nrow = 2)
#' theta <- list(
#'   mu = cbind(c(0, 0), c(1, 0)),
#'   sigma = array(c(a, b), dim = c(2, 2, 2))
#' )
#' compute_average_riemannian(theta)
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
      r <- chol(sigma[,, j])
      # Eigenvalues of R^{-T} Sigma_ell R^{-1}, the congruence of Sigma_ell by
      # Sigma_j^{-1/2}; log-eigenvalues give the AIRM.
      tmp <- backsolve(r, sigma[,, ell], transpose = TRUE)
      ev <- eigen(
        backsolve(r, t(tmp), transpose = TRUE),
        symmetric = TRUE,
        only.values = TRUE
      )$values
      acc <- acc + sqrt(sum(log(pmax(ev, .Machine$double.xmin))^2))
      n_pairs <- n_pairs + 1L
    }
  }
  if (n_pairs < 1L) {
    return(NA_real_)
  }
  acc / n_pairs
}

#' Average pairwise overlap of a Gaussian mixture
#'
#' @description
#' Returns MixSim's **BarOmega**: the unweighted mean of pairwise
#' overlaps
#' \deqn{
#'   \overline{\omega}
#'   =
#'   \frac{2}{J(J-1)}
#'   \sum_{1\le j<\ell\le J}
#'   \bigl(\Omega_{j\ell}+\Omega_{\ell j}\bigr)
#'   \in[0,1]
#' }
#' (up to the MixSim numerical convention), where
#' \eqn{\Omega_{j\ell}=\Pr_{X\sim f_j}(X\text{ classified as }\ell)}
#' already uses the mixture weights \eqn{\boldsymbol{p}} inside the
#' Bayes / MAP rule of [MixSim::overlap()]. Do **not** multiply the
#' directional masses by \eqn{p_j} again.
#'
#' For \eqn{G\ge 4} MixSim's Davies quadrature becomes expensive.
#' The helper then switches to [overlap_gaussian_mc()]: stratified
#' Sobol draws, inverse-transform sampling through precomputed
#' Cholesky factors, and log-density MAP comparisons (`n_mc`
#' draws per component; default 10,000).
#'
#' @param true_theta List with the generative parameters (see `.parse_true_theta()`): `p` (length
#'   \eqn{J} or \eqn{J\times N}), `mu` (\eqn{G\times J}), `sigma`
#'   (\eqn{G\times G\times J}).
#' @param J Number of cell types (components). Defaults to the third
#'   dimension of `sigma`.
#' @param n_mc Monte Carlo draws per component when \eqn{G\ge 4}.
#' @param seed Optional seed forwarded to [overlap_gaussian_mc()].
#' @param verbose If `TRUE` (default), announce the MixSim-to-MC switch
#'   with `cli` when it is installed.
#' @return Scalar average pairwise overlap (MixSim `BarOmega`).
#' @export
#' @seealso [overlap_gaussian_mc()], [compute_average_jeffreys()]
#' @examples
#' set.seed(1)
#' theta <- list(
#'   p = c(0.5, 0.5),
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_overlap(theta)
compute_average_overlap <- function(
  true_theta,
  J = NULL,
  n_mc = 10000L,
  seed = NULL,
  verbose = TRUE
) {
  .overlap_summary(
    true_theta,
    J = J,
    n_mc = n_mc,
    seed = seed,
    verbose = verbose
  )$BarOmega
}

# Mean (`BarOmega`) and maximal (`MaxOmega`) pairwise overlap of a Gaussian
# mixture, from one MixSim call (G < 4) or one Monte Carlo pass (G >= 4).
.overlap_summary <- function(
  true_theta,
  J = NULL,
  n_mc = 10000L,
  seed = NULL,
  verbose = TRUE
) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = TRUE,
    J = J,
    second_moment = "sigma"
  )
  if (theta$G >= 4L) {
    if (isTRUE(verbose)) {
      cli::cli_alert_info(paste(
        "G = {.val {theta$G}} >= 4: MixSim Davies overlap is replaced",
        "by stratified Sobol Monte Carlo ({.fn overlap_gaussian_mc},",
        "n_mc = {.val {n_mc}} draws per component)."
      ))
    }
    return(
      overlap_gaussian_mc(
        true_theta = theta,
        n_mc = n_mc,
        seed = seed,
        J = theta$J
      )
    )
  }
  .check_suggested_package("MixSim", "compute_average_overlap")
  MixSim::overlap(
    Pi = theta$p,
    Mu = t(theta$mu),
    S = theta$sigma
  )
}

#' Average pairwise Jeffreys divergence of a Gaussian mixture
#'
#' @description
#' Returns the \eqn{p_j p_{\ell}}-weighted mean of pairwise Jeffreys
#' (symmetrised KL) divergences between purified Gaussians, with the closed
#' form of the multivariate-normal Kullback--Leibler divergence
#' (<https://statproofbook.github.io/P/mvn-kl.html>)
#' \deqn{
#'   \overline{J}
#'   =
#'   \frac{
#'     \sum_{1\le j<\ell\le J} p_j p_{\ell}\, J(f_{j},f_{\ell})
#'   }{
#'     \sum_{1\le j<\ell\le J} p_j p_{\ell}
#'   }
#'   \in[0,\infty),
#' }
#' with
#' \eqn{f_{j}=\mathcal{N}_{G}(
#'   \boldsymbol{\mu}_{\cdot j},\boldsymbol{\Sigma}_{j}
#' )}.
#' If `p` is omitted it defaults to the equi-balanced vector \eqn{1/J}, which
#' recovers the uniform pairwise average
#' \eqn{2/(J(J-1))\sum_{j<\ell}J(f_{j},f_{\ell})}.
#'
#' @param true_theta List with the generative parameters (see `.parse_true_theta()`): `mu`
#'   (\eqn{G\times J}), `sigma` (\eqn{G\times G\times J}), and optionally `p`
#'   (length \eqn{J} or \eqn{J\times N}). If `p` is missing it is set to
#'   \eqn{(1/J,\ldots,1/J)}.
#' @param J Number of cell types. Defaults to the third dimension of `sigma`.
#'
#' @return Scalar average pairwise Jeffreys divergence.
#' @export
#' @seealso [compute_average_overlap()], [compute_average_riemannian()]
#' @examples
#' set.seed(1)
#' theta <- list(
#'   mu = cbind(c(0, 0), c(3, 0)),
#'   sigma = array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' )
#' compute_average_jeffreys(theta)
compute_average_jeffreys <- function(
  true_theta,
  J = NULL
) {
  theta <- .parse_true_theta(
    true_theta,
    require_p = FALSE,
    J = J,
    second_moment = "sigma"
  )
  p <- theta$p
  if (is.null(p)) {
    p <- rep(1 / theta$J, theta$J)
  }
  mu_gj <- theta$mu
  sigma <- theta$sigma
  J <- theta$J

  pair_sum <- 0
  weight_sum <- 0
  for (j in seq_len(J - 1L)) {
    for (ell in seq.int(j + 1L, J)) {
      w <- p[[j]] * p[[ell]]
      # Jeffreys = KL(f_j || f_l) + KL(f_l || f_j) for two Gaussians.
      theta_j <- solve(sigma[,, j])
      theta_l <- solve(sigma[,, ell])
      delta <- mu_gj[, j] - mu_gj[, ell]
      jeffreys <- 0.5 *
        sum(diag(
          theta_l %*%
            sigma[,, j] +
            theta_j %*% sigma[,, ell] -
            2 * diag(nrow(theta_j))
        )) +
        0.5 * as.numeric(crossprod(delta, (theta_j + theta_l) %*% delta))
      pair_sum <- pair_sum + w * jeffreys
      weight_sum <- weight_sum + w
    }
  }
  if (weight_sum <= 0) {
    stop("Pair weights from `p` are degenerate.")
  }
  pair_sum / weight_sum
}
