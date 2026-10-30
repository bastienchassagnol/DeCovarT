# Interval estimation for interior compositions. Wald (delta method) is
# vcov_ilr_delta() / confint.decovart_fit(). The default of
# confint_decovart() maps a Wald interval for the ILR coordinate back
# through the inverse chart. Profile intervals invert the Wilks
# chi-square_1 calibration. Boundary nulls are in 03_02_boundary_tests.R
# and resampling in 03_03_resampling_inference.R.
# See vignette("01-2-population-inference").

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
        lower = .profile_root(deviance, estimate = p_hat, face = 0),
        upper = .profile_root(deviance, estimate = p_hat, face = 1)
      )
    },
    numeric(3)
  )
  t(intervals)
}

#' Bracket-and-solve one profile-likelihood endpoint
#'
#' Searches between the estimate, where the shifted deviance is
#' \eqn{-\chi^2_{1,1-\alpha}<0}, and a face of the simplex (`0` or `1`).
#' When the deviance stays below the cut-off up to the face, the interval
#' is one-sided and the face is returned. A non-negative deviance at the
#' estimate means the restricted refit failed there, and `NA` is returned
#' with a warning.
#'
#' @noRd
.profile_root <- function(deviance, estimate, face) {
  if (!is.finite(estimate)) {
    return(NA_real_)
  }
  if (abs(face - estimate) < 1e-10) {
    return(face)
  }
  value <- function(x) tryCatch(deviance(x), error = function(e) NA_real_)
  at_estimate <- value(estimate)
  at_face <- value(face)
  if (anyNA(c(at_estimate, at_face))) {
    return(NA_real_)
  }
  if (at_estimate >= 0) {
    warning(
      "The restricted refit at the estimate did not reach the maximum; ",
      "the profile endpoint is set to NA.",
      call. = FALSE
    )
    return(NA_real_)
  }
  if (at_face < 0) {
    return(face)
  }
  tryCatch(
    stats::uniroot(
      deviance,
      lower = min(estimate, face),
      upper = max(estimate, face),
      tol = 1e-8
    )$root,
    error = function(e) NA_real_
  )
}

# Shared interval helpers --------------------------------------------------

#' @noRd
.check_level <- function(level) {
  if (length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1) {
    stop("`level` must lie strictly between 0 and 1.", call. = FALSE)
  }
  level
}

#' @noRd
.level_bounds <- function(level) {
  alpha <- (1 - level) / 2
  paste0(format(100 * c(alpha, 1 - alpha), trim = TRUE), "%")
}

#' @noRd
.percentile_interval <- function(estimates, level) {
  alpha <- (1 - level) / 2
  interval <- t(apply(
    estimates,
    1L,
    stats::quantile,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE,
    names = FALSE
  ))
  colnames(interval) <- .level_bounds(level)
  rownames(interval) <- rownames(estimates)
  interval
}

#' @noRd
.wald_bounds <- function(estimate, se, level, df = Inf) {
  # df = Inf is the normal quantile; a finite df gives a Student-t one.
  z <- stats::qt(1 - (1 - level) / 2, df = df)
  interval <- cbind(estimate - z * se, estimate + z * se)
  colnames(interval) <- .level_bounds(level)
  if (!is.null(names(estimate))) {
    rownames(interval) <- names(estimate)
  }
  interval
}

#' @noRd
.as_bulk_columns <- function(y, mean_signature_matrix) {
  y <- as.matrix(y)
  if (nrow(y) != nrow(mean_signature_matrix)) {
    stop(
      "`y` must have ",
      nrow(mean_signature_matrix),
      " genes (rows), got ",
      nrow(y),
      ".",
      call. = FALSE
    )
  }
  y
}

#' @noRd
.pooled_population_loglik <- function(p, y, mean_signature_matrix, Sigma) {
  sum(vapply(
    seq_len(ncol(y)),
    function(i) {
      loglik_multivariate(p, y[, i], mean_signature_matrix, Sigma)
    },
    numeric(1)
  ))
}

#' @noRd
.pooled_population_score <- function(p, y, mean_signature_matrix, Sigma) {
  scores <- vapply(
    seq_len(ncol(y)),
    function(i) {
      gradient_loglik_unconstrained(p, y[, i], mean_signature_matrix, Sigma)
    },
    numeric(length(p))
  )
  rowSums(matrix(scores, nrow = length(p)))
}

#' @noRd
.pooled_population_hessian <- function(p, y, mean_signature_matrix, Sigma) {
  Reduce(
    `+`,
    lapply(seq_len(ncol(y)), function(i) {
      hessian_loglik_unconstrained(p, y[, i], mean_signature_matrix, Sigma)
    })
  )
}

#' Expected information of the ILR coordinate, for one or more replicates
#'
#' `n_obs` multiplies the per-column Fisher information. Returns `NULL`
#' when `p` is on a face or the ILR information is singular.
#'
#' @noRd
.population_chart_covariance <- function(
  p,
  mean_signature_matrix,
  Sigma,
  n_obs = 1L
) {
  n_obs <- as.integer(n_obs)
  if (length(n_obs) != 1L || is.na(n_obs) || n_obs < 1L) {
    stop("`n_obs` must be a positive integer.", call. = FALSE)
  }
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    return(NULL)
  }
  rho <- isometric_log_ratio(p)
  jac <- jacobian_isometric_logistic(rho)
  info_p <- expected_fisher_unconstrained(p, mean_signature_matrix, Sigma)
  info_rho <- n_obs * (crossprod(jac, info_p %*% jac))
  vcov_rho <- tryCatch(solve(info_rho), error = function(e) NULL)
  if (is.null(vcov_rho)) {
    return(NULL)
  }
  vcov_rho <- .symmetrise(vcov_rho)
  vcov_p <- jac %*% vcov_rho %*% t(jac)
  dimnames(vcov_p) <- list(names(p), names(p))
  list(
    rho = rho,
    jacobian = jac,
    information_rho = .symmetrise(info_rho),
    vcov_rho = vcov_rho,
    vcov_p = vcov_p
  )
}

#' @noRd
.with_seed <- function(seed, expr) {
  if (is.null(seed)) {
    return(expr)
  }
  if (length(seed) != 1L || !is.finite(seed)) {
    stop("`seed` must be a single finite number, or NULL.", call. = FALSE)
  }
  # The RNG kinds are left as they are; the previous state is restored.
  withr::with_seed(seed, expr)
}

#' @noRd
.central_gradient <- function(fn, z, step = 1e-5) {
  h <- step * pmax(1, abs(z))
  centre <- fn(z)
  g <- numeric(length(z))
  for (i in seq_along(z)) {
    forward <- z
    back <- z
    forward[[i]] <- z[[i]] + h[[i]]
    back[[i]] <- z[[i]] - h[[i]]
    upper <- fn(forward)
    lower <- fn(back)
    if (is.finite(upper) && is.finite(lower)) {
      g[[i]] <- (upper - lower) / (2 * h[[i]])
    } else if (is.finite(upper) && is.finite(centre)) {
      g[[i]] <- (upper - centre) / h[[i]]
    } else if (is.finite(lower) && is.finite(centre)) {
      g[[i]] <- (centre - lower) / h[[i]]
    } else {
      g[[i]] <- NA_real_
    }
  }
  g
}

#' @noRd
.empty_proportion_interval <- function(p, level) {
  matrix(
    NA_real_,
    nrow = length(p),
    ncol = 2L,
    dimnames = list(names(p), .level_bounds(level))
  )
}

# ILR back-transformation --------------------------------------------------

#' Push a Wald normal for rho through the inverse ILR chart
#'
#' @noRd
.map_wald_rho <- function(
  rho,
  vcov_rho,
  p_names,
  level,
  n_draws,
  seed,
  df = Inf
) {
  se_rho <- sqrt(pmax(diag(vcov_rho), 0))
  names(se_rho) <- paste0("rho", seq_along(se_rho))
  rho <- stats::setNames(as.numeric(rho), names(se_rho))
  rho_interval <- .wald_bounds(rho, se_rho, level, df = df)
  if (length(rho) == 1L) {
    mapped <- rbind(
      isometric_logistic(rho_interval[[1L]]),
      isometric_logistic(rho_interval[[2L]])
    )
    interval <- cbind(
      pmin(mapped[1L, ], mapped[2L, ]),
      pmax(mapped[1L, ], mapped[2L, ])
    )
    mapping <- "endpoint"
    n_used <- NA_integer_
  } else {
    n_draws <- as.integer(n_draws)
    if (length(n_draws) != 1L || is.na(n_draws) || n_draws < 2L) {
      stop("`n_draws` must be an integer of at least 2.", call. = FALSE)
    }
    # A finite df draws from the multivariate t: normal deviations
    # divided by sqrt(chi^2_df / df).
    draws <- .with_seed(seed, {
      normal <- MASS::mvrnorm(
        n = n_draws,
        mu = rep(0, length(rho)),
        Sigma = vcov_rho
      )
      normal <- matrix(normal, nrow = n_draws)
      if (is.finite(df)) {
        normal <- normal / sqrt(stats::rchisq(n_draws, df = df) / df)
      }
      sweep(normal, 2L, as.numeric(rho), "+")
    })
    p_draws <- t(apply(draws, 1L, isometric_logistic))
    alpha <- (1 - level) / 2
    interval <- t(apply(
      p_draws,
      2L,
      stats::quantile,
      probs = c(alpha, 1 - alpha),
      names = FALSE
    ))
    mapping <- "monte_carlo"
    n_used <- n_draws
  }
  colnames(interval) <- .level_bounds(level)
  rownames(interval) <- p_names
  list(
    rho = rho,
    rho_interval = rho_interval,
    interval = interval,
    mapping = mapping,
    n_draws = n_used
  )
}

#' Confidence intervals by mapping a Wald region through the ILR chart
#'
#' @description
#' Forms a Wald interval for the unconstrained ILR coordinate
#' \eqn{\boldsymbol{\rho}} from the expected Fisher information, then
#' maps it to the simplex by the inverse chart
#' \eqn{\boldsymbol{\psi}}. For two cell types the map is monotone, so
#' the endpoints of the interval for \eqn{\rho} are transformed. For
#' more than two cell types the joint normal draw
#' \eqn{\boldsymbol{\rho}^{(b)}} is transformed and the interval is the
#' marginal quantile of \eqn{\boldsymbol{\psi}(\boldsymbol{\rho}^{(b)})}.
#' Every transformed draw is a composition, so the interval lies in
#' \eqn{(0,1)}.
#'
#' This is the default of [confint_decovart()]. The symmetric delta-method
#' interval on \eqn{\boldsymbol{p}} remains [confint.decovart_fit()].
#' Columns of `y` are replicates that share one composition: their number
#' multiplies the Fisher information.
#'
#' @inheritParams restricted_mle_decovart
#' @param p Composition at which to evaluate the information. `NULL`
#'   uses the pooled maximum likelihood estimate.
#' @param level Confidence level.
#' @param n_draws Monte Carlo size used when there are three or more
#'   cell types. Ignored for two cell types.
#' @param seed Optional seed for the Monte Carlo map. The previous
#'   random-number state is restored. `NULL` leaves the state unchanged
#'   and draws from the current stream.
#'
#' @return A list with `estimate`, `interval` on the simplex, `rho`,
#'   `rho_interval` (component-wise Wald bounds for \eqn{\boldsymbol{\rho}}),
#'   `vcov_rho`, `vcov`, `mapping` (`"endpoint"` or `"monte_carlo"`) and
#'   `n_obs`. On a simplex face the interval is `NA`.
#'
#' @seealso [confint_decovart()], [vcov_ilr_delta()],
#'   [confint.decovart_fit()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' confint_ilr_backtransform(y, mu, Sigma)$interval
confint_ilr_backtransform <- function(
  y,
  mean_signature_matrix,
  Sigma,
  p = NULL,
  level = 0.95,
  n_draws = 2000L,
  seed = NULL,
  epsilon = 10^-8,
  itmax = 500L
) {
  level <- .check_level(level)
  y <- .as_bulk_columns(y, mean_signature_matrix)
  n_obs <- ncol(y)
  if (is.null(p)) {
    p <- restricted_mle_decovart(
      y,
      mean_signature_matrix,
      Sigma,
      epsilon = epsilon,
      itmax = itmax
    )$coefficients
  }
  p <- stats::setNames(as.numeric(p), colnames(mean_signature_matrix))
  chart <- .population_chart_covariance(
    p,
    mean_signature_matrix,
    Sigma,
    n_obs = n_obs
  )
  if (is.null(chart)) {
    warning(
      "ILR back-transformation is undefined on a simplex face.",
      call. = FALSE
    )
    return(list(
      method = "ilr",
      estimate = p,
      interval = .empty_proportion_interval(p, level),
      rho = NA_real_,
      rho_interval = NULL,
      vcov_rho = NULL,
      vcov = matrix(
        NA_real_,
        length(p),
        length(p),
        dimnames = list(names(p), names(p))
      ),
      mapping = NA_character_,
      n_obs = n_obs
    ))
  }
  mapped <- .map_wald_rho(
    rho = chart$rho,
    vcov_rho = chart$vcov_rho,
    p_names = names(p),
    level = level,
    n_draws = n_draws,
    seed = seed
  )
  list(
    method = "ilr",
    estimate = p,
    interval = mapped$interval,
    rho = mapped$rho,
    rho_interval = mapped$rho_interval,
    vcov_rho = chart$vcov_rho,
    vcov = chart$vcov_p,
    mapping = mapped$mapping,
    n_obs = n_obs,
    n_draws = mapped$n_draws
  )
}

# Firth penalisation -------------------------------------------------------

#' @noRd
.jeffreys_penalty_population <- function(z, mean_signature_matrix, Sigma) {
  chart <- .population_chart_covariance(
    isometric_logistic(z),
    mean_signature_matrix,
    Sigma,
    n_obs = 1L
  )
  if (is.null(chart)) {
    return(-Inf)
  }
  ev <- eigen(
    chart$information_rho,
    symmetric = TRUE,
    only.values = TRUE
  )$values
  if (any(!is.finite(ev)) || min(ev) <= 1e-10) {
    return(-Inf)
  }
  0.5 * sum(log(ev))
}

#' Firth-penalised DeCovarT fit
#'
#' @description
#' Maximises
#' \eqn{\ell(\boldsymbol{\rho})+\tfrac12\log\det I_{\boldsymbol{\rho}}}
#' on the ILR chart. The penalty is Jeffreys' prior for the expected
#' Fisher information of one bulk column; extra replicate columns add
#' only a constant, so they do not change the penalty gradient. The
#' returned Wald covariance is that information multiplied by the number
#' of columns and evaluated at the penalised estimate.
#'
#' The penalty is differentiated numerically. The log-likelihood score
#' stays analytic. This is not an exact removal of the \eqn{O(N^{-1})}
#' bias: the convolution is a curved exponential family.
#'
#' @inheritParams confint_ilr_backtransform
#' @param p_start Start on the simplex. `NULL` uses the unpenalised MLE,
#'   pulled into the open simplex when that MLE lies on a face.
#'
#' @return A list with `coefficients`, `rho`, `loglik` (unpenalised),
#'   `penalty`, `penalised_loglik`, `vcov`, `vcov_rho` and `convergence`.
#'
#' @references
#' \insertRef{firthBiasReductionMaximum1993}{DeCovarT}
#'
#' @seealso [confint_decovart()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' fit_firth_decovart(y, mu, Sigma)$coefficients
fit_firth_decovart <- function(
  y,
  mean_signature_matrix,
  Sigma,
  p_start = NULL,
  epsilon = 1e-8,
  itmax = 200L
) {
  y <- .as_bulk_columns(y, mean_signature_matrix)
  nms <- colnames(mean_signature_matrix)
  if (is.null(p_start)) {
    p_start <- restricted_mle_decovart(
      y,
      mean_signature_matrix,
      Sigma,
      epsilon = epsilon,
      itmax = max(as.integer(itmax), 50L)
    )$coefficients
  }
  p_start <- repair_simplex(p_start, open = TRUE, nms = nms)
  z0 <- isometric_log_ratio(p_start)
  penalty <- function(z) {
    .jeffreys_penalty_population(z, mean_signature_matrix, Sigma)
  }
  objective <- function(z) {
    p_try <- isometric_logistic(z)
    ll <- .pooled_population_loglik(p_try, y, mean_signature_matrix, Sigma)
    pen <- penalty(z)
    if (!is.finite(ll) || !is.finite(pen)) {
      return(-1e10)
    }
    ll + pen
  }
  gradient <- function(z) {
    p_try <- isometric_logistic(z)
    g_ll <- .ilr_gradient(
      z,
      .pooled_population_score(p_try, y, mean_signature_matrix, Sigma)
    )
    g_pen <- .central_gradient(penalty, z)
    if (any(!is.finite(g_pen))) {
      g_pen <- rep(0, length(z))
    }
    g_ll + g_pen
  }
  fit <- stats::optim(
    par = z0,
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(
      fnscale = -1,
      reltol = epsilon,
      maxit = as.integer(itmax)
    )
  )
  p_hat <- repair_simplex(isometric_logistic(fit$par), open = TRUE, nms = nms)
  chart <- .population_chart_covariance(
    p_hat,
    mean_signature_matrix,
    Sigma,
    n_obs = ncol(y)
  )
  pen <- penalty(fit$par)
  list(
    coefficients = p_hat,
    rho = stats::setNames(fit$par, paste0("rho", seq_along(fit$par))),
    loglik = .pooled_population_loglik(
      p_hat,
      y,
      mean_signature_matrix,
      Sigma
    ),
    penalty = pen,
    penalised_loglik = fit$value,
    vcov = if (is.null(chart)) {
      matrix(
        NA_real_,
        length(p_hat),
        length(p_hat),
        dimnames = list(nms, nms)
      )
    } else {
      chart$vcov_p
    },
    vcov_rho = if (is.null(chart)) NULL else chart$vcov_rho,
    convergence = fit$convergence
  )
}

# Dispatcher ---------------------------------------------------------------

#' Choose a confidence interval for a pooled DeCovarT composition
#'
#' @description
#' One entry point for the interior intervals of the population
#' convolution (covariance weights \eqn{p_j^2}, no technical scales).
#' The default, `method = "ilr"`, is [confint_ilr_backtransform()].
#'
#' \describe{
#'   \item{`ilr`}{Wald interval for \eqn{\boldsymbol{\rho}}, mapped by
#'     the inverse ILR chart.}
#'   \item{`wald`}{Delta method on \eqn{\boldsymbol{p}}. Endpoints can
#'     leave \eqn{[0,1]}.}
#'   \item{`firth`}{Wald interval at the Jeffreys-penalised estimate
#'     ([fit_firth_decovart()]). `interval_ilr` is the same estimate
#'     mapped through the inverse chart.}
#'   \item{`profile`}{Likelihood-ratio interval
#'     ([confint_profile_decovart()]). This is the interior calibration.
#'     A null on a face is [lrt_decovart()], not an interval.}
#'   \item{`bootstrap`}{Parametric percentile interval
#'     ([bootstrap_decovart()]).}
#'   \item{`reference`}{Resample the reference, then refit
#'     ([reference_bootstrap_decovart()]). Two or more donors in every
#'     cell type use the donor bootstrap and keep the bulk fixed.
#'     Otherwise columns are resampled within each cell type and the
#'     bulk is regenerated. Cell-type labels are not shuffled; that
#'     check is [equivariance_check_decovart()].}
#' }
#'
#' Columns of `y` share one composition. [fit_decovart()] instead
#' estimates a separate composition per column, and
#' [confint.decovart_fit()] is the per-column delta method.
#'
#' @inheritParams confint_ilr_backtransform
#' @param method Interval to compute. See Description.
#' @param n_boot Bootstrap size for `method = "bootstrap"` or
#'   `"reference"`.
#' @param reference_profiles,donor_ids Passed to
#'   [reference_bootstrap_decovart()] when `method = "reference"`.
#' @param resample `"donors"` or `"cells"`. `NULL` chooses from
#'   `donor_ids` as in Description.
#' @param regenerate_bulk Whether the reference bootstrap redraws the
#'   bulk. `NULL` follows the automatic rule in Description.
#' @param ... Passed to the selected method (`epsilon`, `itmax`,
#'   `celltypes`, `null_value`, `ridge`, `n_replicates`).
#'
#' @return A list with `method`, `estimate` and `interval`. Other
#'   elements depend on `method`.
#'
#' @seealso [lrt_decovart()], [confint_single_cell()]
#' @family decovart_inference
#' @export
#' @examples
#' mu <- matrix(c(20, 40, 15, 40, 20, 25), nrow = 3)
#' colnames(mu) <- paste0("ct", 1:2)
#' Sigma <- array(c(diag(3), diag(3)), dim = c(3, 3, 2))
#' y <- drop(mu %*% c(0.6, 0.4))
#' confint_decovart(y, mu, Sigma)$interval
#' confint_decovart(y, mu, Sigma, method = "wald")$interval
confint_decovart <- function(
  y,
  mean_signature_matrix = NULL,
  Sigma = NULL,
  method = c(
    "ilr",
    "wald",
    "firth",
    "profile",
    "bootstrap",
    "reference"
  ),
  level = 0.95,
  p = NULL,
  n_draws = 2000L,
  n_boot = 500L,
  seed = NULL,
  reference_profiles = NULL,
  donor_ids = NULL,
  resample = NULL,
  regenerate_bulk = NULL,
  ...
) {
  method <- match.arg(method)
  level <- .check_level(level)
  if (!identical(method, "reference")) {
    if (is.null(mean_signature_matrix) || is.null(Sigma)) {
      stop(
        "`mean_signature_matrix` and `Sigma` are required unless ",
        "method = \"reference\".",
        call. = FALSE
      )
    }
  }
  if (identical(method, "ilr")) {
    return(confint_ilr_backtransform(
      y,
      mean_signature_matrix,
      Sigma,
      p = p,
      level = level,
      n_draws = n_draws,
      seed = seed,
      ...
    ))
  }
  if (identical(method, "wald")) {
    y_mat <- .as_bulk_columns(y, mean_signature_matrix)
    if (is.null(p)) {
      p <- restricted_mle_decovart(
        y_mat,
        mean_signature_matrix,
        Sigma,
        ...
      )$coefficients
    }
    p <- stats::setNames(as.numeric(p), colnames(mean_signature_matrix))
    chart <- .population_chart_covariance(
      p,
      mean_signature_matrix,
      Sigma,
      n_obs = ncol(y_mat)
    )
    if (is.null(chart)) {
      warning(
        "Wald interval is undefined on a simplex face.",
        call. = FALSE
      )
      interval <- .empty_proportion_interval(p, level)
    } else {
      se <- sqrt(pmax(diag(chart$vcov_p), 0))
      names(se) <- names(p)
      interval <- .wald_bounds(p, se, level)
    }
    return(list(
      method = "wald",
      estimate = p,
      interval = interval,
      vcov = if (is.null(chart)) NULL else chart$vcov_p,
      n_obs = ncol(y_mat)
    ))
  }
  if (identical(method, "firth")) {
    fit <- fit_firth_decovart(
      y,
      mean_signature_matrix,
      Sigma,
      p_start = p,
      ...
    )
    se <- sqrt(pmax(diag(fit$vcov), 0))
    names(se) <- names(fit$coefficients)
    ilr <- confint_ilr_backtransform(
      y,
      mean_signature_matrix,
      Sigma,
      p = fit$coefficients,
      level = level,
      n_draws = n_draws,
      seed = seed
    )
    return(list(
      method = "firth",
      estimate = fit$coefficients,
      interval = .wald_bounds(fit$coefficients, se, level),
      interval_ilr = ilr$interval,
      fit = fit
    ))
  }
  if (identical(method, "profile")) {
    profile <- confint_profile_decovart(
      y,
      mean_signature_matrix,
      Sigma,
      level = level,
      ...
    )
    interval <- profile[, c("lower", "upper"), drop = FALSE]
    colnames(interval) <- .level_bounds(level)
    return(list(
      method = "profile",
      estimate = stats::setNames(profile[, "estimate"], rownames(profile)),
      interval = interval
    ))
  }
  if (identical(method, "bootstrap")) {
    boot <- .with_seed(
      seed,
      bootstrap_decovart(
        y,
        mean_signature_matrix,
        Sigma,
        p = p,
        n_boot = n_boot,
        level = level,
        ...
      )
    )
    return(list(
      method = "bootstrap",
      estimate = boot$p_simulated,
      interval = boot$interval,
      bootstrap = boot
    ))
  }
  if (is.null(reference_profiles)) {
    stop(
      "method = \"reference\" needs `reference_profiles`.",
      call. = FALSE
    )
  }
  plan <- .reference_resample_plan(donor_ids)
  if (is.null(resample)) {
    resample <- plan$method
  } else {
    resample <- match.arg(resample, c("donors", "cells"))
  }
  if (is.null(regenerate_bulk)) {
    regenerate_bulk <- if (identical(resample, "donors")) {
      FALSE
    } else {
      TRUE
    }
  }
  ref <- .with_seed(
    seed,
    reference_bootstrap_decovart(
      y = y,
      reference_profiles = reference_profiles,
      donor_ids = donor_ids,
      method = resample,
      n_boot = n_boot,
      level = level,
      regenerate_bulk = regenerate_bulk,
      p = p,
      ...
    )
  )
  list(
    method = "reference",
    resample = ref$method,
    regenerate_bulk = isTRUE(regenerate_bulk),
    estimate = ref$p_hat,
    interval = ref$interval,
    reference = ref
  )
}
