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
  global <- tensor::tensor(as.numeric(weights), Sigma, alongA = 1L, alongB = 3L)
  (global + t(global)) / 2
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
