# Covariance factorisation and Newton driver shared by the population
# convolution (weights p_j^2) and the single-cell mixture (weights p_j).

#' Weighted sum of cell-type covariances
#'
#' \eqn{\sum_j w_j\boldsymbol{\Sigma}_j}, symmetrised. The population
#' convolution uses \eqn{w_j=p_j^2}, the single-cell mixture
#' \eqn{w_j=p_j}.
#'
#' @noRd
.mixture_covariance <- function(weights, Sigma) {
  .symmetrise(tensor::tensor(
    as.numeric(weights),
    Sigma,
    alongA = 1L,
    alongB = 3L
  ))
}

#' Symmetric part of a square matrix
#'
#' @noRd
.symmetrise <- function(x) {
  x <- as.matrix(x)
  (x + t(x)) / 2
}

#' Cached Cholesky factorisation of a weighted covariance mixture
#'
#' Newton-type solvers evaluate the log-likelihood, score and Hessian at
#' the same trial point within one iteration. This single-slot cache,
#' keyed on exact equality of `weights` and `Sigma`, factorises the
#' \eqn{G\times G} mixture once and returns the matrix, its upper
#' Cholesky factor, the log-determinant and the precision.
#'
#' @noRd
.mixture_factorisation <- local({
  cache <- new.env(parent = emptyenv())
  cache$value <- NULL

  function(weights, Sigma) {
    cached <- cache$value
    if (
      !is.null(cached) &&
        identical(cached$weights, weights) &&
        identical(cached$Sigma, Sigma)
    ) {
      return(cached$value)
    }
    mixture <- .ensure_spd(.mixture_covariance(weights, Sigma))
    chol_factor <- base::chol(mixture)
    value <- list(
      matrix = mixture,
      chol = chol_factor,
      log_det = 2 * sum(log(diag(chol_factor))),
      inverse = chol2inv(chol_factor)
    )
    cache$value <- list(weights = weights, Sigma = Sigma, value = value)
    value
  }
})

#' Trust-region Newton ascent in ILR coordinates
#'
#' Maximises `objective` with [stats::nlminb()], which, given an analytic
#' Hessian, runs the PORT model/trust-region Newton routine: each step
#' maximises the quadratic model inside a ball whose radius adapts to the
#' agreement between predicted and actual increase. There is no line
#' search and no fixed step length.
#'
#' `abs.tol = 0` switches off PORT's absolute function test
#' \eqn{|f|<\mathrm{abs.tol}}, which is meaningful only for a nonnegative
#' objective; a negative log-likelihood can cross zero.
#'
#' @param z0 Start in ILR coordinates.
#' @param objective,gradient,hessian Functions of `z` returning the
#'   log-likelihood, its score and its Hessian (to be maximised).
#' @param epsilon Relative function and step tolerances.
#' @param itmax Maximum number of Newton iterations.
#'
#' @return List with `par`, `loglik`, `convergence` (0 means one of the
#'   PORT convergence tests was met), `iterations` and `message`.
#'
#' @noRd
.newton_ilr <- function(z0, objective, gradient, hessian, epsilon, itmax) {
  fit <- stats::nlminb(
    start = z0,
    objective = function(z) -objective(z),
    gradient = function(z) -gradient(z),
    hessian = function(z) -hessian(z),
    control = list(
      iter.max = itmax,
      rel.tol = epsilon,
      x.tol = epsilon,
      xf.tol = epsilon,
      abs.tol = 0
    )
  )
  message <- if (is.null(fit$message)) "" else as.character(fit$message)
  if (!is.finite(fit$objective)) {
    warning(
      "Newton-Raphson did not return a finite objective; the start is kept.",
      call. = FALSE
    )
    return(list(
      par = z0,
      loglik = objective(z0),
      convergence = 1L,
      iterations = 0L,
      message = "non-finite objective"
    ))
  }
  list(
    par = fit$par,
    loglik = -fit$objective,
    convergence = as.integer(fit$convergence),
    iterations = as.integer(fit$iterations),
    message = message
  )
}

#' Maximum likelihood on a face of the simplex
#'
#' Coordinates in `fixed` are substituted, so an exact zero is
#' representable; the free block is \eqn{(1-\sum\mathrm{fixed})
#' \boldsymbol{\psi}(\boldsymbol{z})} on the ILR chart of the free face.
#' That block is fitted by [.newton_ilr()] from its barycentre, with the
#' score and Hessian in the full composition pulled back by the chain
#' rule. The trust region keeps the first step bounded: a quasi-Newton
#' step scaled by a large pooled score can otherwise jump to a vertex,
#' where the likelihood is flat, and stop there.
#'
#' @param nms Cell-type names, length \eqn{J}.
#' @param fixed Named ratios held fixed, or `NULL`.
#' @param loglik,score,hessian Functions of the full composition
#'   \eqn{\boldsymbol{p}} (length \eqn{J}) returning the log-likelihood,
#'   its gradient and its Hessian.
#' @param epsilon,itmax Passed to [.newton_ilr()].
#'
#' @return A list with `coefficients`, `loglik`, `fixed` and
#'   `convergence` (`code`, `iterations`, `message`, `rdm`).
#'
#' @noRd
.fit_on_face <- function(nms, fixed, loglik, score, hessian, epsilon, itmax) {
  n_celltypes <- length(nms)
  idx <- match(names(fixed), nms)
  if (length(fixed) > 0L && (is.null(names(fixed)) || anyNA(idx))) {
    stop(
      "`fixed` must be named by columns of `mean_signature_matrix`.",
      call. = FALSE
    )
  }
  fixed_value <- as.numeric(fixed)
  if (
    length(fixed_value) > 0L &&
      (any(!is.finite(fixed_value)) ||
        any(fixed_value < 0) ||
        sum(fixed_value) > 1 + 1e-12)
  ) {
    stop(
      "`fixed` must hold non-negative ratios summing to at most 1.",
      call. = FALSE
    )
  }
  free <- setdiff(seq_len(n_celltypes), idx)
  remaining <- max(0, 1 - sum(fixed_value))
  assemble <- function(z) {
    p <- numeric(n_celltypes)
    p[idx] <- fixed_value
    if (length(free) == 1L) {
      p[free] <- remaining
    } else if (length(free) > 1L && remaining > 0) {
      p[free] <- remaining * isometric_logistic(z)
    }
    stats::setNames(p, nms)
  }
  result <- function(p, value, convergence) {
    list(
      coefficients = p,
      loglik = value,
      fixed = stats::setNames(fixed_value, nms[idx]),
      convergence = convergence
    )
  }
  if (length(free) <= 1L || remaining <= 0) {
    p_hat <- assemble(numeric(0))
    return(result(
      p_hat,
      loglik(p_hat),
      list(
        code = 0L,
        iterations = 0L,
        message = "closed form",
        rdm = NA_real_
      )
    ))
  }
  fit <- .newton_ilr(
    rep(0, length(free) - 1L),
    objective = function(z) loglik(assemble(z)),
    gradient = function(z) {
      .ilr_gradient(z, remaining * score(assemble(z))[free])
    },
    hessian = function(z) {
      p <- assemble(z)
      .ilr_hessian(
        z,
        remaining * score(p)[free],
        remaining^2 * hessian(p)[free, free, drop = FALSE]
      )
    },
    epsilon = epsilon,
    itmax = itmax
  )
  result(
    assemble(fit$par),
    fit$loglik,
    list(
      code = fit$convergence,
      iterations = fit$iterations,
      message = fit$message,
      rdm = NA_real_
    )
  )
}
