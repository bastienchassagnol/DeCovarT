# Single-cell linear mixture: Sigma(p) = sum_j p_j Sigma_j, with two
# technical scales c (mean) and phi (covariance), profiled out in closed
# form. The ILR chart (R/02_00_ilr_chart.R) and the covariance cache and
# trust-region Newton driver (R/02_00_likelihood_core.R) are shared with
# the population convolution (R/02_02_population_decovart_estimate_ratios.R),
# whose covariance uses p_j^2 and is not called here.

#' @noRd
.check_single_cell_inputs <- function(y, mean_signature_matrix, Sigma) {
  y <- as.numeric(y)
  mean_signature_matrix <- as.matrix(mean_signature_matrix)
  n_genes <- nrow(mean_signature_matrix)
  n_types <- ncol(mean_signature_matrix)
  if (n_genes <= n_types) {
    stop(
      "The two technical scales need more genes than cell types ",
      "(G > J); got G = ",
      n_genes,
      " and J = ",
      n_types,
      ".",
      call. = FALSE
    )
  }
  if (length(y) != n_genes || any(!is.finite(y))) {
    stop("`y` must have one finite entry per gene.", call. = FALSE)
  }
  if (
    length(dim(Sigma)) != 3L ||
      !identical(dim(Sigma)[1:2], c(n_genes, n_genes)) ||
      dim(Sigma)[[3L]] != n_types
  ) {
    stop(
      "`Sigma` must be a G x G x J array matching the mean.",
      call. = FALSE
    )
  }
  if (any(!is.finite(mean_signature_matrix)) || any(!is.finite(Sigma))) {
    stop("The reference moments must be finite.", call. = FALSE)
  }
  list(y = y, mean = mean_signature_matrix, n_types = n_types)
}

#' @noRd
.positive_scalar <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || !is.finite(x) || x <= 0) {
    stop("`", name, "` must be a positive finite scalar.", call. = FALSE)
  }
  x
}

#' Free-scale nonnegative least squares
#'
#' Solves \eqn{\min_{w\ge 0}\lVert y-\mu w\rVert_2^2}, then
#' \eqn{c=\mathbf{1}^{\top}w} and \eqn{p=w/c}. This is the quadratic
#' programme of DeconRNASeq with a free scale. Supports are enumerated,
#' which is exact for the modest \eqn{J} of a cell-type panel.
#'
#' @noRd
.nnls_free_scale <- function(mean_signature_matrix, y, tol = 1e-8) {
  y <- as.numeric(y)
  n_types <- ncol(mean_signature_matrix)
  if (n_types > 12L) {
    stop(
      "The free-scale start enumerates supports and is limited ",
      "to 12 cell types.",
      call. = FALSE
    )
  }
  best <- rep(0, n_types)
  best_sse <- sum(y * y)
  n_masks <- 2^n_types
  for (mask in seq_len(n_masks - 1L)) {
    active <- as.logical(intToBits(mask)[seq_len(n_types)])
    coef <- stats::lm.fit(
      mean_signature_matrix[, active, drop = FALSE],
      y
    )$coefficients
    if (anyNA(coef) || any(coef < -tol)) {
      next
    }
    coef[coef < 0] <- 0
    weights <- rep(0, n_types)
    weights[active] <- coef
    residual <- y - drop(mean_signature_matrix %*% weights)
    sse <- sum(residual * residual)
    if (sse < best_sse) {
      best_sse <- sse
      best <- weights
    }
  }
  total <- sum(best)
  if (!is.finite(total) || total <= 0) {
    p <- rep(1 / n_types, n_types)
    total <- NA_real_
  } else {
    p <- best / total
  }
  names(p) <- colnames(mean_signature_matrix)
  list(weights = best, ols_scale = total, p = p)
}

#' Closed-form scales at a fixed composition
#'
#' For fixed \eqn{p}, the log-likelihood is maximised by the generalised
#' least-squares mean scale \eqn{\hat c=(m^\top Wm)^{-1}m^\top Wy} and the
#' Mahalanobis variance scale \eqn{\hat\phi=r^\top Wr/G}, with
#' \eqn{m=\tilde\mu p}, \eqn{W=S_W(p)^{-1}} and \eqn{r=y-\hat c m}. These
#' are the two moment conditions of the article and the profile maximisers.
#' With `strict = FALSE`, a non-positive scale returns `NULL` instead of an
#' error, so a trial point of the optimiser can be rejected.
#'
#' @noRd
.single_cell_scales_at_p <- function(
  y,
  mean_signature_matrix,
  Sigma,
  p,
  strict = TRUE
) {
  fail <- function(msg) {
    if (strict) {
      stop(msg, call. = FALSE)
    }
    NULL
  }
  weight <- .mixture_factorisation(as.numeric(p), Sigma)$inverse
  mean_vector <- drop(mean_signature_matrix %*% p)
  weighted_mean <- drop(weight %*% mean_vector)
  denominator <- sum(mean_vector * weighted_mean)
  if (!is.finite(denominator) || denominator <= 0) {
    return(fail("The GLS denominator for the mean scale is not positive."))
  }
  mean_scale <- sum(weighted_mean * y) / denominator
  if (!is.finite(mean_scale) || mean_scale <= 0) {
    return(fail("The GLS mean scale is not positive."))
  }
  residual <- y - mean_scale * mean_vector
  variance_scale <- sum(residual * drop(weight %*% residual)) / length(y)
  if (!is.finite(variance_scale) || variance_scale <= 0) {
    return(fail("The Mahalanobis variance scale is not positive."))
  }
  list(mean_scale = mean_scale, variance_scale = variance_scale)
}

#' Log-likelihood of one single-cell bulk sample
#'
#' @description
#' Gaussian log-likelihood of
#' \eqn{\boldsymbol{y}\sim\mathcal{N}(c\tilde{\boldsymbol{\mu}}\boldsymbol{p},
#' \phi\boldsymbol{S}_W(\boldsymbol{p}))}
#' with
#' \eqn{\boldsymbol{S}_W(\boldsymbol{p})=\sum_j p_j\boldsymbol{\Sigma}_{W,j}}.
#' The scales \eqn{c} and \eqn{\phi} are held fixed. Up to
#' \eqn{-\tfrac{G}{2}\log(2\pi)},
#' \deqn{
#'   \ell
#'   =
#'   -\tfrac{G}{2}\log\phi
#'   -\tfrac{1}{2}\log\det\boldsymbol{S}_W(\boldsymbol{p})
#'   -\tfrac{1}{2\phi}
#'   \boldsymbol{r}^{\top}\boldsymbol{S}_W(\boldsymbol{p})^{-1}\boldsymbol{r},
#' }
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-c\tilde{\boldsymbol{\mu}}\boldsymbol{p}}.
#' The population objective [loglik_multivariate()] is a different model:
#' its covariance is \eqn{\sum_j p_j^2\boldsymbol{\Sigma}_j}.
#'
#' @param p Composition of length \eqn{J}. It does not have to sum to 1
#'   (numerical derivatives step off the simplex).
#' @param y Bulk vector of length \eqn{G}.
#' @param mean_signature_matrix Reference mean, \eqn{G\times J}.
#' @param Sigma Cell-level covariances, \eqn{G\times G\times J}.
#' @param mean_scale Positive plug-in \eqn{c=an}.
#' @param variance_scale Positive plug-in \eqn{\phi=a^2 n}.
#'
#' @return Scalar log-likelihood.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' p <- c(0.7, 0.3)
#' y <- 20 * drop(mu %*% p)
#' loglik_single_cell(p, y, mu, sigma, 20, 4)
#' @export
loglik_single_cell <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  y <- as.numeric(y)
  p <- as.numeric(p)
  mean_scale <- .positive_scalar(mean_scale, "mean_scale")
  variance_scale <- .positive_scalar(variance_scale, "variance_scale")
  factor <- .mixture_factorisation(as.numeric(p), Sigma)
  residual <- y - mean_scale * drop(mean_signature_matrix %*% p)
  whitened <- backsolve(factor$chol, residual, transpose = TRUE)
  -0.5 *
    length(y) *
    log(variance_scale) -
    0.5 * factor$log_det -
    0.5 / variance_scale * sum(whitened * whitened)
}

#' Score of the single-cell log-likelihood in \eqn{\boldsymbol{p}}
#'
#' @description
#' With \eqn{\boldsymbol{W}=\boldsymbol{S}_W(\boldsymbol{p})^{-1}} and
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-c\tilde{\boldsymbol{\mu}}\boldsymbol{p}},
#' \deqn{
#'   U_j
#'   =
#'   \phi^{-1}\bigl(
#'     c\,\tilde{\boldsymbol{\mu}}_{\cdot j}^{\top}\boldsymbol{W}\boldsymbol{r}
#'     +\tfrac{1}{2}
#'     \boldsymbol{r}^{\top}\boldsymbol{W}\boldsymbol{\Sigma}_j
#'     \boldsymbol{W}\boldsymbol{r}
#'   \bigr)
#'   -\tfrac{1}{2}\mathrm{tr}(\boldsymbol{W}\boldsymbol{\Sigma}_j).
#' }
#' This is the score in the article. It is not
#' [gradient_loglik_unconstrained()], whose trace term carries an extra
#' \eqn{p_j} because \eqn{\partial\boldsymbol{\Sigma}/\partial p_j
#' =2p_j\boldsymbol{\Sigma}_j}.
#'
#' @inheritParams loglik_single_cell
#' @return Numeric vector of length \eqn{J}.
#' @keywords internal
#' @export
gradient_loglik_single_cell <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  y <- as.numeric(y)
  p <- as.numeric(p)
  mean_scale <- .positive_scalar(mean_scale, "mean_scale")
  variance_scale <- .positive_scalar(variance_scale, "variance_scale")
  factor <- .mixture_factorisation(as.numeric(p), Sigma)
  weight <- factor$inverse
  residual <- y - mean_scale * drop(mean_signature_matrix %*% p)
  whitened <- drop(weight %*% residual)
  n_types <- length(p)
  trace_term <- vapply(
    seq_len(n_types),
    function(j) sum(weight * Sigma[,, j]),
    numeric(1)
  )
  quad_term <- vapply(
    seq_len(n_types),
    function(j) sum(whitened * drop(Sigma[,, j] %*% whitened)),
    numeric(1)
  )
  mean_term <- drop(crossprod(mean_signature_matrix, whitened))
  as.numeric(
    -0.5 *
      trace_term +
      (mean_scale / variance_scale) * mean_term +
      (0.5 / variance_scale) * quad_term
  )
}

#' Observed Hessian of the single-cell log-likelihood in \eqn{\boldsymbol{p}}
#'
#' @description
#' Second derivatives of [loglik_single_cell()] with respect to
#' \eqn{\boldsymbol{p}}. Newton-Raphson uses this observed Hessian.
#' The expected information in the article,
#' \eqn{I_{jk}=c^2\phi^{-1}
#' \tilde{\boldsymbol{\mu}}_{\cdot j}^{\top}\boldsymbol{W}
#' \tilde{\boldsymbol{\mu}}_{\cdot k}
#' +\tfrac{1}{2}\mathrm{tr}(
#' \boldsymbol{W}\boldsymbol{\Sigma}_j\boldsymbol{W}\boldsymbol{\Sigma}_k)},
#' is \eqn{\mathbb{E}[-\mathbf{H}]} and is not what [numDeriv::hessian()]
#' computes, so it is not the matrix passed to the solver.
#'
#' @inheritParams loglik_single_cell
#' @return Symmetric numeric \eqn{J\times J} matrix.
#' @keywords internal
#' @export
hessian_loglik_single_cell <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  y <- as.numeric(y)
  p <- as.numeric(p)
  mean_scale <- .positive_scalar(mean_scale, "mean_scale")
  variance_scale <- .positive_scalar(variance_scale, "variance_scale")
  factor <- .mixture_factorisation(as.numeric(p), Sigma)
  weight <- factor$inverse
  residual <- y - mean_scale * drop(mean_signature_matrix %*% p)
  whitened <- drop(weight %*% residual)
  n_types <- length(p)
  n_genes <- length(y)
  sigma_w <- matrix(0, n_genes, n_types)
  for (j in seq_len(n_types)) {
    sigma_w[, j] <- drop(Sigma[,, j] %*% whitened)
  }
  weighted_sigma_w <- weight %*% sigma_w
  weighted_mean <- weight %*% mean_signature_matrix
  weighted_sigma <- lapply(
    seq_len(n_types),
    function(j) weight %*% Sigma[,, j]
  )
  trace_prod <- matrix(0, n_types, n_types)
  for (j in seq_len(n_types)) {
    for (k in j:n_types) {
      tr <- sum(weighted_sigma[[j]] * t(weighted_sigma[[k]]))
      trace_prod[j, k] <- tr
      trace_prod[k, j] <- tr
    }
  }
  mean_block <- crossprod(mean_signature_matrix, weighted_mean)
  quad <- crossprod(sigma_w, weighted_sigma_w)
  cross_mean <- crossprod(mean_signature_matrix, weighted_sigma_w)
  cross_sigma <- crossprod(sigma_w, weighted_mean)
  ratio <- mean_scale / variance_scale
  hessian <- 0.5 *
    trace_prod -
    (mean_scale * ratio) * mean_block -
    quad / variance_scale -
    ratio * cross_mean -
    ratio * cross_sigma
  unname(.symmetrise(hessian))
}

#' Single-cell log-likelihood in ILR coordinates
#'
#' @inheritParams loglik_single_cell
#' @param z ILR coordinate of length \eqn{J-1}.
#' @param V Optional ILR basis. See [isometric_logistic()].
#' @return Scalar log-likelihood.
#' @keywords internal
#' @export
loglik_single_cell_ilr <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  loglik_single_cell(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    mean_scale,
    variance_scale
  )
}

#' Score of the single-cell log-likelihood in ILR coordinates
#'
#' Chain rule \eqn{\mathbf{J}_{\psi}^{\top}U} through
#' [jacobian_isometric_logistic()].
#'
#' @inheritParams loglik_single_cell_ilr
#' @return Numeric vector of length \eqn{J-1}.
#' @keywords internal
#' @export
gradient_loglik_single_cell_ilr <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  .ilr_gradient(
    z,
    gradient_loglik_single_cell(
      p,
      y,
      mean_signature_matrix,
      Sigma,
      mean_scale,
      variance_scale
    ),
    V
  )
}

#' Observed Hessian of the single-cell log-likelihood in ILR coordinates
#'
#' Second-order chain rule through [jacobian_isometric_logistic()] and
#' [hessian_isometric_logistic()], as for [hessian_loglik_constrained()].
#'
#' @inheritParams loglik_single_cell_ilr
#' @return Symmetric numeric \eqn{(J-1)\times(J-1)} matrix.
#' @keywords internal
#' @export
hessian_loglik_single_cell_ilr <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  args <- list(p, y, mean_signature_matrix, Sigma, mean_scale, variance_scale)
  .ilr_hessian(
    z,
    do.call(gradient_loglik_single_cell, args),
    do.call(hessian_loglik_single_cell, args),
    V
  )
}

#' Stage-i scales for one single-cell bulk sample
#'
#' @description
#' Generalised method of moments for the two technical scales, at a
#' composition that is either supplied or taken from the free-scale
#' DeconRNASeq start
#' \eqn{(c^{(0)},\boldsymbol{p}^{(0)})
#' =\arg\min_{c>0,\,p\in\Delta}\lVert y-c\tilde\mu p\rVert_2^2}.
#' The start is floored into the open simplex. At that composition,
#' \eqn{\hat c
#' =(\boldsymbol{m}^{\top}\boldsymbol{W}\boldsymbol{m})^{-1}
#' \boldsymbol{m}^{\top}\boldsymbol{W}\boldsymbol{y}}
#' and
#' \eqn{\hat\phi=\boldsymbol{r}^{\top}\boldsymbol{W}\boldsymbol{r}/G},
#' with \eqn{\boldsymbol{W}=\boldsymbol{S}_W(\boldsymbol{p})^{-1}}.
#' Then \eqn{\hat a=\hat\phi/\hat c} and \eqn{\hat n=\hat c^2/\hat\phi}.
#' \eqn{\hat n} is a cell count only when the two platforms share their
#' scale; it is not rounded.
#'
#' \eqn{\hat c} is not the library-size ratio
#' \eqn{\mathbf{1}^{\top}y/\sum_j p_j S_j}.
#'
#' @inheritParams loglik_single_cell
#' @param p Optional composition. `NULL` runs the free-scale start.
#'
#' @return List with `p`, `ols_scale` (the least-squares \eqn{c^{(0)}},
#'   or `NA` when `p` was supplied), `mean_scale`, `variance_scale`,
#'   `platform_scale` and `n_cells`.
#'
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' set.seed(1)
#' y <- 30 * drop(mu %*% c(0.6, 0.4)) + rnorm(4, sd = 0.5)
#' single_cell_nuisance_scales(y, mu, sigma)
#' @export
#' @seealso [deconvolute_single_cell()],
#'   [deconvolute_ratios_deconrnaseq()]
single_cell_nuisance_scales <- function(
  y,
  mean_signature_matrix,
  Sigma,
  p = NULL
) {
  checked <- .check_single_cell_inputs(y, mean_signature_matrix, Sigma)
  y <- checked$y
  mean_signature_matrix <- checked$mean
  nms <- colnames(mean_signature_matrix)
  ols_scale <- NA_real_
  if (is.null(p)) {
    start <- .nnls_free_scale(mean_signature_matrix, y)
    p <- start$p
    ols_scale <- start$ols_scale
  }
  p <- repair_simplex(p, open = TRUE, nms = nms)
  scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p)
  mean_scale <- scales$mean_scale
  variance_scale <- scales$variance_scale
  list(
    p = p,
    ols_scale = ols_scale,
    mean_scale = mean_scale,
    variance_scale = variance_scale,
    platform_scale = variance_scale / mean_scale,
    n_cells = mean_scale^2 / variance_scale
  )
}

#' Observed Hessian blocks in \eqn{(\boldsymbol{p},c,\phi)}
#'
#' `hessian_pp` is the composition block with the scales held fixed.
#' `h_peta` is the \eqn{J\times 2} cross block and `h_etaeta` the
#' \eqn{2\times 2} block of \eqn{\eta=(c,\phi)}. The profile Hessian in
#' \eqn{\boldsymbol{p}} is the Schur complement of `h_etaeta`.
#'
#' @noRd
.single_cell_hessian_blocks <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  hessian_pp <- hessian_loglik_single_cell(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    mean_scale,
    variance_scale
  )
  weight <- .mixture_factorisation(as.numeric(p), Sigma)$inverse
  mean_vector <- drop(mean_signature_matrix %*% p)
  residual <- y - mean_scale * mean_vector
  whitened <- drop(weight %*% residual)
  weighted_mean <- drop(weight %*% mean_vector)
  sigma_w <- vapply(
    seq_along(p),
    function(j) drop(Sigma[,, j] %*% whitened),
    numeric(length(y))
  )
  sigma_w <- matrix(sigma_w, nrow = length(y))
  mean_w <- drop(crossprod(mean_signature_matrix, whitened))
  quad <- colSums(sigma_w * whitened)
  phi <- variance_scale
  mean_wm <- drop(crossprod(mean_signature_matrix, weighted_mean))
  sigma_wm <- drop(crossprod(sigma_w, weighted_mean))
  h_pc <- (mean_w - mean_scale * mean_wm - sigma_wm) / phi
  h_pphi <- -(mean_scale * mean_w + 0.5 * quad) / phi^2
  h_cc <- -sum(mean_vector * weighted_mean) / phi
  h_cphi <- -sum(mean_vector * whitened) / phi^2
  h_phiphi <- length(y) / (2 * phi^2) - sum(residual * whitened) / phi^3
  list(
    hessian_pp = hessian_pp,
    h_peta = cbind(h_pc, h_pphi),
    h_etaeta = matrix(c(h_cc, h_cphi, h_cphi, h_phiphi), 2L, 2L)
  )
}

#' Profile Hessian of the single-cell log-likelihood in \eqn{\boldsymbol{p}}
#'
#' Exact second derivative of
#' \eqn{\ell^{\ast}(p)=\ell(p,\hat c(p),\hat\phi(p))}: the Schur complement
#' \eqn{H_{pp}-H_{p\eta}H_{\eta\eta}^{-1}H_{\eta p}} of the observed
#' Hessian in \eqn{(p,\eta)}, \eqn{\eta=(c,\phi)}, evaluated at the
#' closed-form scales. With \eqn{w=Wr},
#' \eqn{H_{p_jc}=\phi^{-1}(\tilde\mu_j^\top w-c\tilde\mu_j^\top Wm
#' -m^\top W\Sigma_j w)},
#' \eqn{H_{p_j\phi}=-\phi^{-2}(c\tilde\mu_j^\top w+\tfrac12 w^\top\Sigma_j w)},
#' \eqn{H_{cc}=-m^\top Wm/\phi}, \eqn{H_{c\phi}=-m^\top w/\phi^2} and
#' \eqn{H_{\phi\phi}=G/(2\phi^2)-r^\top w/\phi^3}.
#'
#' @noRd
.profile_hessian_single_cell <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale,
  variance_scale
) {
  blocks <- .single_cell_hessian_blocks(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    mean_scale,
    variance_scale
  )
  out <- blocks$hessian_pp -
    blocks$h_peta %*% solve(blocks$h_etaeta, t(blocks$h_peta))
  unname(.symmetrise(out))
}

#' Profile log-likelihood of one single-cell bulk sample in ILR coordinates
#'
#' @description
#' \eqn{\ell^{\ast}(\boldsymbol{\rho})
#' =\ell(\boldsymbol{\psi}(\boldsymbol{\rho}),\hat c,\hat\phi)}, with the
#' scales recomputed in closed form at every composition
#' (see [single_cell_nuisance_scales()]). By the envelope theorem its score
#' is the score of [loglik_single_cell()] at those scales, mapped by
#' \eqn{\mathbf{J}_{\psi}^{\top}}. Its Hessian is the Schur complement of
#' the observed Hessian in \eqn{(\boldsymbol{p},c,\phi)}, mapped by the
#' second-order chain rule.
#'
#' The profile objective returns `-Inf` where a closed-form scale is not
#' positive, so the trust-region driver rejects that trial point.
#'
#' @inheritParams loglik_single_cell_ilr
#' @return `profile_loglik_single_cell_ilr()`: scalar.
#'   `gradient_profile_loglik_single_cell_ilr()`: vector of length
#'   \eqn{J-1}. `hessian_profile_loglik_single_cell_ilr()`: symmetric
#'   \eqn{(J-1)\times(J-1)} matrix.
#' @keywords internal
#' @export
profile_loglik_single_cell_ilr <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  scales <- .single_cell_scales_at_p(
    y,
    mean_signature_matrix,
    Sigma,
    p,
    strict = FALSE
  )
  if (is.null(scales)) {
    return(-Inf)
  }
  loglik_single_cell(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    scales$mean_scale,
    scales$variance_scale
  )
}

#' @rdname profile_loglik_single_cell_ilr
#' @export
gradient_profile_loglik_single_cell_ilr <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p)
  .ilr_gradient(
    z,
    gradient_loglik_single_cell(
      p,
      y,
      mean_signature_matrix,
      Sigma,
      scales$mean_scale,
      scales$variance_scale
    ),
    V
  )
}

#' @rdname profile_loglik_single_cell_ilr
#' @export
hessian_profile_loglik_single_cell_ilr <- function(
  z,
  y,
  mean_signature_matrix,
  Sigma,
  V = NULL
) {
  p <- isometric_logistic(z, V)
  scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p)
  args <- list(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    scales$mean_scale,
    scales$variance_scale
  )
  .ilr_hessian(
    z,
    do.call(gradient_loglik_single_cell, args),
    do.call(.profile_hessian_single_cell, args),
    V
  )
}

#' Observed information of one single-cell fit
#'
#' @description
#' At a fitted composition, the covariance of the ILR coordinate
#' \eqn{\hat{\boldsymbol{\rho}}} is the inverse of the negative profile
#' Hessian [hessian_profile_loglik_single_cell_ilr()]. The delta method
#' carries it to \eqn{\hat{\boldsymbol{p}}} with
#' \eqn{\mathbf{J}_{\psi}}. The two technical scales use the other Schur
#' complement of the joint observed Hessian in
#' \eqn{(\boldsymbol{\rho},c,\phi)}: estimating the composition widens
#' their variances, and \eqn{I_{c\phi}} is not assumed to make that
#' complement diagonal.
#'
#' Both matrices are observed information, not expected Fisher
#' information. They treat the reference moments as known.
#'
#' @inheritParams loglik_single_cell
#' @param p Fitted composition. A boundary coordinate makes the ILR
#'   chart undefined.
#' @param mean_scale,variance_scale Scales at `p`. `NULL` recomputes the
#'   closed forms.
#' @param scales_known If `TRUE`, the scales are treated as known: the
#'   composition covariance inverts the fixed-scale observed Hessian
#'   [hessian_loglik_single_cell_ilr()] at `mean_scale` and
#'   `variance_scale`, and `vcov_nuisance` is `NA`. This is the
#'   covariance of a `"fixed_scales"` fit. For a `"two_stage"` fit it
#'   ignores the uncertainty of the start scales, so its intervals are
#'   too narrow.
#' @param warn If `TRUE`, warn when either information matrix is not
#'   positive definite.
#'
#' @return A list with `rho`, `vcov_rho`, `vcov` (delta method for
#'   \eqn{\boldsymbol{p}}), `mean_scale`, `variance_scale` and
#'   `vcov_nuisance` (order \eqn{c}, \eqn{\phi}). Singular blocks are
#'   filled with `NA`.
#'
#' @seealso [deconvolute_single_cell()], [confint_single_cell()]
#' @export
#' @examples
#' mu <- matrix(c(2, 0.2, 0.2, 2, 0.2, 3, 0.2, 3), nrow = 4)
#' sigma <- array(rep(diag(0.05, 4), 2), c(4, 4, 2))
#' y <- 20 * drop(mu %*% c(0.7, 0.3)) + c(0.3, -0.2, 0.1, -0.1)
#' vcov_single_cell(c(0.7, 0.3), y, mu, sigma)
vcov_single_cell <- function(
  p,
  y,
  mean_signature_matrix,
  Sigma,
  mean_scale = NULL,
  variance_scale = NULL,
  scales_known = FALSE,
  warn = TRUE
) {
  checked <- .check_single_cell_inputs(y, mean_signature_matrix, Sigma)
  y <- checked$y
  mean_signature_matrix <- checked$mean
  nms <- colnames(mean_signature_matrix)
  p <- repair_simplex(p, open = TRUE, nms = nms)
  if (is.null(mean_scale) || is.null(variance_scale)) {
    scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p)
    mean_scale <- scales$mean_scale
    variance_scale <- scales$variance_scale
  } else {
    mean_scale <- .positive_scalar(mean_scale, "mean_scale")
    variance_scale <- .positive_scalar(variance_scale, "variance_scale")
  }
  rho <- isometric_log_ratio(p)
  jac <- jacobian_isometric_logistic(rho)
  na_rho <- matrix(NA_real_, length(rho), length(rho))
  na_p <- matrix(NA_real_, length(p), length(p), dimnames = list(nms, nms))
  na_eta <- matrix(
    NA_real_,
    2L,
    2L,
    dimnames = list(
      c("mean_scale", "variance_scale"),
      c("mean_scale", "variance_scale")
    )
  )
  h_rho <- .symmetrise(hessian_loglik_single_cell_ilr(
    rho,
    y,
    mean_signature_matrix,
    Sigma,
    mean_scale,
    variance_scale
  ))
  info_rho <- if (isTRUE(scales_known)) {
    -h_rho
  } else {
    .symmetrise(
      -hessian_profile_loglik_single_cell_ilr(
        rho,
        y,
        mean_signature_matrix,
        Sigma
      )
    )
  }
  vcov_rho <- .invert_pd(
    info_rho,
    warn = warn,
    what = if (isTRUE(scales_known)) "fixed-scale" else "profile"
  )
  if (is.null(vcov_rho)) {
    vcov_rho <- na_rho
    vcov_p <- na_p
  } else {
    vcov_p <- jac %*% vcov_rho %*% t(jac)
    dimnames(vcov_p) <- list(nms, nms)
  }
  blocks <- .single_cell_hessian_blocks(
    p,
    y,
    mean_signature_matrix,
    Sigma,
    mean_scale,
    variance_scale
  )
  h_cross <- crossprod(jac, blocks$h_peta)
  vcov_eta <- tryCatch(
    {
      schur <- .symmetrise(
        blocks$h_etaeta - t(h_cross) %*% solve(h_rho, h_cross)
      )
      solved <- .invert_pd(-schur, warn = warn, what = "nuisance")
      if (is.null(solved)) na_eta else solved
    },
    error = function(e) {
      if (isTRUE(warn)) {
        warning(
          "Nuisance information could not be inverted.",
          call. = FALSE
        )
      }
      na_eta
    }
  )
  if (isTRUE(scales_known)) {
    vcov_eta <- na_eta
  }
  dimnames(vcov_eta) <- dimnames(na_eta)
  list(
    rho = rho,
    vcov_rho = vcov_rho,
    vcov = vcov_p,
    mean_scale = mean_scale,
    variance_scale = variance_scale,
    vcov_nuisance = vcov_eta
  )
}


#' @noRd
.invert_pd <- function(info, warn, what) {
  ev <- eigen(info, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(ev)) || min(ev) <= 1e-8) {
    if (isTRUE(warn)) {
      warning(
        "The ",
        what,
        " information is not positive definite.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  .symmetrise(solve(info))
}

#' Single-cell deconvolution of one bulk sample
#'
#' @description
#' Fits one bulk sample from the linear single-cell mixture
#' \eqn{\boldsymbol{y}\sim\mathcal{N}(c\tilde{\boldsymbol{\mu}}\boldsymbol{p},
#' \phi\boldsymbol{S}_W(\boldsymbol{p}))}.
#'
#' The start is [single_cell_nuisance_scales()]: nonnegative least squares
#' with a free scale (DeconRNASeq's quadratic programme without the simplex
#' equality), then the closed-form scales at that composition.
#'
#' With `profile = TRUE` (default) the two scales are profiled out: for
#' every trial composition they are recomputed in closed form, so the
#' weight matrix \eqn{\boldsymbol{S}_W(\boldsymbol{p})^{-1}} is rebuilt at
#' every iterate, and Newton-Raphson runs on the profile log-likelihood of
#' the \eqn{J-1} ILR coordinates ([profile_loglik_single_cell_ilr()]). Its
#' limit is the joint maximum-likelihood estimate of
#' \eqn{(\boldsymbol{p},c,\phi)}, reached in one loop.
#'
#' With `profile = FALSE` the scales stay at their start values (a
#' two-stage plug-in fit), which is not the joint maximum. Supplying both
#' `mean_scale` and `variance_scale` fixes them at those values.
#'
#' Every fit uses the analytic score and observed Hessian with the PORT
#' trust region of [stats::nlminb()]: the step maximises the quadratic
#' model inside a ball whose radius grows or shrinks with the ratio of
#' actual to predicted increase; there is no line search and no fixed
#' step length. The returned `vcov` is the delta-method image of the
#' inverse negative profile Hessian; intervals are
#' [confint_single_cell()].
#'
#' @inheritParams loglik_single_cell
#' @param epsilon Relative function and step tolerances passed to
#'   [stats::nlminb()].
#' @param itmax Maximum number of Newton iterations.
#' @param profile If `TRUE`, profile the two scales out; if `FALSE`, keep
#'   them at their start values.
#' @param initial_p Optional open-simplex start. `NULL` uses the free-scale
#'   least-squares start.
#' @param mean_scale,variance_scale Optional fixed scales. Supply both
#'   or neither.
#' @param estimand `cellular` (default) interprets the fitted
#'   composition as a cell fraction, which is valid when the mean
#'   keeps the cell sizes \eqn{S_j}. `transcript` interprets it as a
#'   transcript fraction; pass `cell_size` to convert it back to a
#'   cell fraction.
#' @param cell_size Optional \eqn{S_j}. Defaults to the column sums of
#'   the mean, which equal the cell sizes for a `cellular_raw`
#'   reference.
#'
#' @return A list with `composition` (the fitted simplex coordinate),
#'   `p`, `q`, `rho`, `estimand`, `fit_type` (`"profile"`, `"two_stage"`
#'   or `"fixed_scales"`), `mean_scale`, `variance_scale`,
#'   `platform_scale`, `n_cells`, `cell_size`, `ols_scale`, `loglik`,
#'   `convergence` (0 from [stats::nlminb()] means success), `iterations`,
#'   `p_start`, `message`, `vcov_rho`, `vcov` and `vcov_nuisance`
#'   (see [vcov_single_cell()]). A singular information matrix is
#'   returned as `NA` without failing the fit.
#'
#' @examples
#' genes <- paste0("g", 1:4)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(
#'   c(2, 2, 0.2, 0.2, 0.2, 0.2, 3, 3),
#'   nrow = 4,
#'   dimnames = list(genes, cts)
#' )
#' sigma <- array(0, c(4, 4, 2), dimnames = list(genes, genes, cts))
#' sigma[,, 1] <- diag(0.05, 4)
#' sigma[,, 2] <- diag(0.05, 4)
#' set.seed(1)
#' sim <- simulate_bulk_from_cells(
#'   mu,
#'   sigma,
#'   p = c(0.8, 0.2),
#'   n_cells = 2000,
#'   platform_scale = 1
#' )
#' fit <- deconvolute_single_cell(sim$Y[, 1], mu, sigma)
#' round(fit$p, 3)
#' @export
#' @seealso [single_cell_nuisance_scales()],
#'   [profile_loglik_single_cell_ilr()],
#'   [deconvolute_ratios_Newton_Raphson()],
#'   [rescale_single_cell_moments()]
deconvolute_single_cell <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 1e-8,
  itmax = 100L,
  profile = TRUE,
  initial_p = NULL,
  mean_scale = NULL,
  variance_scale = NULL,
  estimand = c("cellular", "transcript"),
  cell_size = NULL
) {
  estimand <- match.arg(estimand)
  checked <- .check_single_cell_inputs(y, mean_signature_matrix, Sigma)
  y <- checked$y
  mean_signature_matrix <- checked$mean
  nms <- colnames(mean_signature_matrix)
  fixed <- !is.null(mean_scale) || !is.null(variance_scale)
  if (fixed && (is.null(mean_scale) || is.null(variance_scale))) {
    stop(
      "`mean_scale` and `variance_scale` must be supplied together.",
      call. = FALSE
    )
  }
  ols_scale <- NA_real_
  if (is.null(initial_p)) {
    start <- single_cell_nuisance_scales(y, mean_signature_matrix, Sigma)
    p_start <- start$p
    ols_scale <- start$ols_scale
    start_scales <- start[c("mean_scale", "variance_scale")]
  } else {
    p_start <- repair_simplex(initial_p, open = TRUE, nms = nms)
    start_scales <- if (fixed) {
      NULL
    } else {
      .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p_start)
    }
  }
  fit_type <- if (fixed) {
    "fixed_scales"
  } else if (isTRUE(profile)) {
    "profile"
  } else {
    "two_stage"
  }
  if (fixed) {
    mean_scale <- .positive_scalar(mean_scale, "mean_scale")
    variance_scale <- .positive_scalar(variance_scale, "variance_scale")
  } else {
    mean_scale <- start_scales$mean_scale
    variance_scale <- start_scales$variance_scale
  }

  z_start <- isometric_log_ratio(p_start)
  if (identical(fit_type, "profile")) {
    fit <- .newton_ilr(
      z_start,
      objective = function(z) {
        profile_loglik_single_cell_ilr(z, y, mean_signature_matrix, Sigma)
      },
      gradient = function(z) {
        gradient_profile_loglik_single_cell_ilr(
          z,
          y,
          mean_signature_matrix,
          Sigma
        )
      },
      hessian = function(z) {
        hessian_profile_loglik_single_cell_ilr(
          z,
          y,
          mean_signature_matrix,
          Sigma
        )
      },
      epsilon = epsilon,
      itmax = as.integer(itmax)
    )
  } else {
    args <- list(
      y = y,
      mean_signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      mean_scale = mean_scale,
      variance_scale = variance_scale
    )
    fit <- .newton_ilr(
      z_start,
      objective = function(z) do.call(loglik_single_cell_ilr, c(list(z), args)),
      gradient = function(z) {
        do.call(gradient_loglik_single_cell_ilr, c(list(z), args))
      },
      hessian = function(z) {
        do.call(hessian_loglik_single_cell_ilr, c(list(z), args))
      },
      epsilon = epsilon,
      itmax = as.integer(itmax)
    )
  }
  p <- repair_simplex(isometric_logistic(fit$par), open = TRUE, nms = nms)
  if (identical(fit_type, "profile")) {
    scales <- .single_cell_scales_at_p(y, mean_signature_matrix, Sigma, p)
    mean_scale <- scales$mean_scale
    variance_scale <- scales$variance_scale
  }
  fitted_composition <- p

  if (is.null(cell_size)) {
    cell_size <- colSums(mean_signature_matrix)
  } else {
    cell_size <- as.numeric(cell_size)
    if (
      length(cell_size) != length(p) ||
        any(!is.finite(cell_size)) ||
        any(cell_size <= 0)
    ) {
      stop("`cell_size` must be a positive length-J vector.", call. = FALSE)
    }
  }
  names(cell_size) <- nms
  if (estimand == "cellular") {
    q <- p * cell_size
    q <- q / sum(q)
  } else {
    q <- p
    ratio <- q / cell_size
    p <- ratio / sum(ratio)
  }
  names(p) <- nms
  names(q) <- nms
  uncertainty <- tryCatch(
    vcov_single_cell(
      fitted_composition,
      y,
      mean_signature_matrix,
      Sigma,
      mean_scale,
      variance_scale,
      scales_known = !identical(fit_type, "profile"),
      warn = FALSE
    ),
    error = function(e) NULL
  )

  list(
    composition = fitted_composition,
    p = p,
    q = q,
    rho = isometric_log_ratio(fitted_composition),
    estimand = estimand,
    fit_type = fit_type,
    mean_scale = mean_scale,
    variance_scale = variance_scale,
    platform_scale = variance_scale / mean_scale,
    n_cells = mean_scale^2 / variance_scale,
    cell_size = cell_size,
    ols_scale = ols_scale,
    loglik = loglik_single_cell(
      fitted_composition,
      y,
      mean_signature_matrix,
      Sigma,
      mean_scale,
      variance_scale
    ),
    convergence = fit$convergence,
    iterations = fit$iterations,
    p_start = p_start,
    message = fit$message,
    vcov_rho = if (is.null(uncertainty)) NULL else uncertainty$vcov_rho,
    vcov = if (is.null(uncertainty)) NULL else uncertainty$vcov,
    vcov_nuisance = if (is.null(uncertainty)) {
      NULL
    } else {
      uncertainty$vcov_nuisance
    }
  )
}
