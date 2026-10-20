# Interval estimation for interior compositions: the Cramer--Rao / ILR
# delta-method (Wald) covariance stored by fit_decovart(), and profile-
# likelihood intervals that invert the Wilks chi-square_1 calibration. Both
# assume the true composition lies in the open simplex; boundary nulls are
# handled in 03_02_boundary_tests.R.
# See vignette("theory-DeCovarT-MLE-properties").

# Wald intervals -----------------------------------------------------------

#' Expected Fisher information of unconstrained \eqn{\boldsymbol{p}}
#'
#' @description
#' For the multivariate-normal mean--covariance map of the DeCovarT
#' convolution,
#' \eqn{\boldsymbol{y}\sim\mathcal{N}_{G}(\boldsymbol{\mu}\boldsymbol{p},
#' \boldsymbol{\Sigma}(\boldsymbol{p}))} with
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_j p_j^{2}\boldsymbol{\Sigma}_j}
#' and precision
#' \eqn{\boldsymbol{\Omega}(\boldsymbol{p})=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}},
#' the expected Fisher information has entries
#' \deqn{
#'   I(\boldsymbol{p})_{jk}
#'   =
#'   \boldsymbol{\mu}_{\cdot j}^{\top}
#'   \boldsymbol{\Omega}(\boldsymbol{p})
#'   \boldsymbol{\mu}_{\cdot k}
#'   +
#'   2 p_j p_k\,
#'   \mathrm{tr}\bigl(
#'     \boldsymbol{\Omega}(\boldsymbol{p})\boldsymbol{\Sigma}_j
#'     \boldsymbol{\Omega}(\boldsymbol{p})\boldsymbol{\Sigma}_k
#'   \bigr).
#' }
#' The first summand is the mean contribution
#' (an \eqn{\boldsymbol{\Omega}}-inner product of signature columns); the
#' second is the covariance contribution of the quadratic map
#' \eqn{\boldsymbol{p}\mapsto\boldsymbol{\Sigma}(\boldsymbol{p})}. See the
#' multivariate-normal formula on
#' <https://en.wikipedia.org/wiki/Fisher_information#Multivariate_normal_distribution>.
#'
#' @param p Numeric proportions on the open simplex.
#' @param mean_signature_matrix Mean signature \eqn{\boldsymbol{\mu}}
#'   (\eqn{G\times J}).
#' @param Sigma Cell-type covariances \eqn{G\times G\times J}.
#'
#' @return Symmetric \eqn{J\times J} expected Fisher information matrix
#'   \eqn{I(\boldsymbol{p})}.
#'
#' @seealso [vcov_ilr_delta()], [vcov.decovart_fit()],
#'   [confint.decovart_fit()], [hessian_loglik_unconstrained()]
#'
#' @export
#' @examples
#' # Two cell types, bivariate Gaussians: I(p) is 2 x 2.
#' p <- c(0.6, 0.4)
#' mu <- cbind(c(0, 0), c(3, 0))
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' I <- expected_fisher_unconstrained(p, mu, Sigma)
#' I
expected_fisher_unconstrained <- function(
  p,
  mean_signature_matrix,
  Sigma
) {
  blocks <- .information_blocks(p, mean_signature_matrix, Sigma)
  blocks$mean + blocks$cov
}

#' Cramer--Rao / ILR delta-method covariance of \eqn{\hat{\boldsymbol{p}}}
#'
#' @description
#' Maps the expected Fisher information of unconstrained proportions
#' through the isometric log-ratio (ILR) chart and back to the simplex
#' via the delta method.
#'
#' Let \eqn{\boldsymbol{p}=\operatorname{softmax}(\mathbf{V}\boldsymbol{z})}
#' with Jacobian
#' \eqn{\mathbf{J}_{\boldsymbol{\psi}}=\mathbf{S}(\boldsymbol{p})
#' \mathbf{V}}
#' ([jacobian_isometric_logistic()]). Fisher information transforms as
#' \deqn{
#'   I_{\boldsymbol{z}}
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}
#'   I(\boldsymbol{p})
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   =
#'   \mathbf{V}^{\mathsf{T}}
#'   \mathbf{S}(\boldsymbol{p})
#'   I(\boldsymbol{p})
#'   \mathbf{S}(\boldsymbol{p})
#'   \mathbf{V}.
#' }
#' The first-order delta method then yields
#' \deqn{
#'   \mathrm{Var}(\hat{\boldsymbol{p}})
#'   \approx
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   I_{\boldsymbol{z}}^{-1}
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}.
#' }
#' This simplex covariance is invariant to orthogonal rotations of
#' \eqn{\mathbf{V}}. The construction is undefined on the simplex
#' boundary (the log-ratio chart blows up); the function then returns
#' `NA` with a warning.
#'
#' Wald intervals in [vcov.decovart_fit()], [confint.decovart_fit()]
#' and [deconvolute_ratios()] all use this expected-Fisher ILR map.
#' The observed Hessian [hessian_loglik_constrained()] is for Newton
#' steps only; it is not the Wald covariance.
#'
#' @inheritParams expected_fisher_unconstrained
#' @param warn If `FALSE`, skip the boundary / singularity warnings
#'   (Monte Carlo loops). Default `TRUE`.
#'
#' @return Symmetric \eqn{J\times J} asymptotic covariance of
#'   \eqn{\hat{\boldsymbol{p}}}, or a matrix of `NA` if the bound is
#'   undefined / singular.
#'
#' @seealso [expected_fisher_unconstrained()], [vcov.decovart_fit()],
#'   [confint.decovart_fit()], [jacobian_isometric_logistic()]
#'
#' @export
#' @examples
#' # Two cell types: a single unconstrained ILR coordinate maps back to
#' # a rank-1 covariance on the simplex (rows sum to zero).
#' p <- c(ct1 = 0.6, ct2 = 0.4)
#' mu <- cbind(ct1 = c(0, 0), ct2 = c(3, 0))
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' v <- vcov_ilr_delta(p, mu, Sigma)
#' round(v, 4)
#' abs(rowSums(v)) < 1e-10
vcov_ilr_delta <- function(
  p,
  mean_signature_matrix,
  Sigma,
  warn = TRUE
) {
  nms <- names(p)
  n_celltypes <- length(p)
  out <- matrix(
    NA_real_,
    n_celltypes,
    n_celltypes,
    dimnames = list(nms, nms)
  )
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    if (isTRUE(warn)) {
      warning(
        "Proportions on the simplex boundary; Wald vcov is undefined.",
        call. = FALSE
      )
    }
    return(out)
  }
  info_p <- expected_fisher_unconstrained(p, mean_signature_matrix, Sigma)
  z <- isometric_log_ratio(p)
  jac <- jacobian_isometric_logistic(z)
  info_z <- t(jac) %*% info_p %*% jac
  vcov_z <- tryCatch(
    solve(info_z),
    error = function(e) {
      if (isTRUE(warn)) {
        warning(
          "Expected Fisher information in ILR coordinates is singular.",
          call. = FALSE
        )
      }
      NULL
    }
  )
  if (is.null(vcov_z)) {
    return(out)
  }
  vcov_p <- jac %*% vcov_z %*% t(jac)
  dimnames(vcov_p) <- list(nms, nms)
  vcov_p
}

#' @noRd
.ilr_wald_se <- function(
  p,
  mean_signature_matrix,
  Sigma,
  warn = FALSE
) {
  nms <- names(p)
  if (is.null(nms)) {
    nms <- colnames(mean_signature_matrix)
  }
  v <- tryCatch(
    vcov_ilr_delta(
      p,
      mean_signature_matrix,
      Sigma,
      warn = warn
    ),
    error = function(e) NULL
  )
  if (is.null(v)) {
    return(stats::setNames(rep(NA_real_, length(p)), nms))
  }
  stats::setNames(sqrt(pmax(diag(as.matrix(v)), 0)), nms)
}

# Profile-likelihood intervals ---------------------------------------------

#' Profile-likelihood confidence intervals for cellular ratios
#'
#' @description
#' Inverts the likelihood-ratio test of [lrt_decovart()] coordinate by
#' coordinate: the interval collects every \eqn{c} with
#' \eqn{2\{\ell(\hat{\boldsymbol{p}})-\ell^{\mathrm{prof}}_j(c)\}
#' \le\chi^{2}_{1,1-\alpha}}. Because likelihood ratios are invariant
#' under reparametrisation, the interval respects
#' \eqn{0\le p_j\le 1} without the delta-method linearisation that makes
#' Wald intervals unreliable near a simplex face.
#'
#' @inheritParams restricted_mle_decovart
#' @param level Confidence level \eqn{1-\alpha}.
#' @param celltypes Names of the cell types to profile (columns of
#'   `mean_signature_matrix`); defaults to all.
#'
#' @return A matrix with one row per cell type and columns `estimate`,
#'   `lower`, `upper`.
#'
#' @seealso [confint.decovart_fit()] for the Wald analogue,
#'   [lrt_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' confint_profile_decovart(y, mu, Sigma)
confint_profile_decovart <- function(
  y,
  mean_signature_matrix,
  Sigma,
  level = 0.95,
  celltypes = colnames(mean_signature_matrix),
  epsilon = 10^-8,
  itmax = 500L
) {
  if (!all(celltypes %in% colnames(mean_signature_matrix))) {
    stop(
      "`celltypes` must be columns of `mean_signature_matrix`.",
      call. = FALSE
    )
  }
  full <- restricted_mle_decovart(
    y,
    mean_signature_matrix,
    Sigma,
    epsilon = epsilon,
    itmax = itmax
  )
  cutoff <- stats::qchisq(level, df = 1L)
  intervals <- vapply(
    celltypes,
    function(celltype) {
      p_hat <- full$coefficients[[celltype]]
      deviance <- function(value) {
        restricted <- restricted_mle_decovart(
          y,
          mean_signature_matrix,
          Sigma,
          fixed = stats::setNames(value, celltype),
          epsilon = epsilon,
          itmax = itmax
        )
        2 * (full$loglik - restricted$loglik) - cutoff
      }
      c(
        estimate = p_hat,
        lower = .profile_root(deviance, lower = 0, upper = p_hat),
        upper = .profile_root(deviance, lower = p_hat, upper = 1)
      )
    },
    numeric(3)
  )
  t(intervals)
}

#' Bracket-and-solve one profile-likelihood endpoint
#'
#' Returns the boundary itself when the deviance never crosses the
#' cut-off inside the bracket (a genuinely one-sided interval).
#'
#' @noRd
.profile_root <- function(deviance, lower, upper) {
  if (!is.finite(lower) || !is.finite(upper) || upper - lower < 1e-10) {
    return(if (upper - lower < 1e-10) lower else NA_real_)
  }
  f_lower <- tryCatch(deviance(lower), error = function(e) NA_real_)
  f_upper <- tryCatch(deviance(upper), error = function(e) NA_real_)
  if (anyNA(c(f_lower, f_upper))) {
    return(NA_real_)
  }
  if (f_lower * f_upper > 0) {
    # Deviance stays below the cut-off across the bracket: the profile
    # interval runs into the simplex face.
    return(if (f_lower < 0) lower else upper)
  }
  tryCatch(
    stats::uniroot(deviance, lower = lower, upper = upper, tol = 1e-8)$root,
    error = function(e) NA_real_
  )
}
