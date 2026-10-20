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
#' @examples
#' sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' scaled <- scale_covariance_array(sigma, 4)
#' scaled[, , 1]
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
#' @keywords internal
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
