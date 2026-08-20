# Algebraic operations ------------------------------------------------------

#' Additive logistic transform (unconstrained coordinates to the simplex)
#'
#' @description
#' Implements the reparametrisation
#' \eqn{\boldsymbol{\psi}:\boldsymbol{\rho}\mapsto\boldsymbol{p}} used in the
#' article, sending unconstrained coordinates
#' \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}} to cellular proportions
#' \eqn{\boldsymbol{p}\in\Delta^{J-1}}. This is the *additive logistic
#' transform* of Aitchison, i.e. the inverse additive log-ratio map
#' (\eqn{\mathrm{alr}^{-1}}), equivalently a softmax with the last category
#' \eqn{J} pinned as reference (\eqn{\rho_J\equiv 0}).
#'
#' @details
#' With \eqn{A=\sum_{k=1}^{J-1}\mathrm{e}^{\rho_k}+1},
#' \deqn{
#'   p_j=\frac{\mathrm{e}^{\rho_j}}{A}\quad(j<J),\qquad
#'   p_J=\frac{1}{A}.
#' }
#' Equivalently,
#' \eqn{\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})} with
#' \eqn{\boldsymbol{\psi}(\boldsymbol{\rho})\propto
#' (\mathrm{e}^{\rho_1},\ldots,\mathrm{e}^{\rho_{J-1}},1)^{\mathsf{T}}}.
#' Jacobians and Hessians of both maps are derived in the package vignette
#' `vignette("softmax-alr-derivatives", package = "DeCovarT")`.
#' See also [compositions::alrInv()].
#'
#' @param rho Numeric vector \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}} of
#'   unconstrained additive log-ratio coordinates (reference cell type
#'   \eqn{J}).
#'
#' @return Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}} on the unit
#'   simplex (\eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1},
#'   \eqn{\boldsymbol{p}\ge\mathbf{0}}).
#'
#' @seealso The inverse map (additive log-ratio) is documented as
#'   `additive_log_ratio()` on this help page.
#'
#' @examples
#' rho <- c(0.2, -0.5)
#' p <- additive_logistic(rho)
#' sum(p)
#' additive_log_ratio(p)
#' @export
additive_logistic <- function(rho) {
  p <- c(exp(rho[seq_along(rho)]), 1)
  return(p / sum(p))
}

#' Additive log-ratio transform \eqn{\boldsymbol{p}\mapsto\boldsymbol{\rho}}
#'
#' Recovers the unconstrained additive log-ratio coordinates
#' \eqn{\rho_j=\ln(p_j/p_J)} for \eqn{j=1,\ldots,J-1}, with the last part
#' \eqn{p_J} as reference. This is Aitchison's additive log-ratio
#' (\eqn{\mathrm{alr}}) transform, equivalently the multinomial-logit link
#' with reference category \eqn{J} (see [compositions::alr()] and
#' `vignette("softmax-alr-derivatives", package = "DeCovarT")`).
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}} on the open simplex.
#'
#' @return Numeric vector \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}}.
#'
#' @rdname additive_logistic
#' @export
additive_log_ratio <- function(p) {
  num_cells <- length(p)
  return(log(p[1:num_cells - 1] / p[num_cells]))
}

#' Evaluate a matrix-induced bilinear form
#'
#' Computes \eqn{\boldsymbol{x}^{\mathsf{T}}\boldsymbol{A}\boldsymbol{y}}.
#' When \eqn{\boldsymbol{A}} is symmetric positive definite (as for a
#' non-degenerate Gaussian covariance or precision), this is the
#' \eqn{\boldsymbol{A}}-inner product. Prefer this name over "dot product",
#' which is reserved for \eqn{\boldsymbol{x}^{\mathsf{T}}\boldsymbol{y}}.
#'
#' Implementation uses [base::crossprod()] as
#' `drop(crossprod(x, A %*% y))`, which is the standard efficient route to
#' a bilinear / quadratic form in R (avoids an explicit transpose of
#' \eqn{\boldsymbol{x}} and a temporary outer product).
#'
#' @param x Numeric vector.
#' @param A Numeric square matrix, compatible with `x` and `y`.
#' @param y Numeric vector of the same length as `x` (default `x`,
#'   which yields the quadratic form
#'   \eqn{\boldsymbol{x}^{\mathsf{T}}\boldsymbol{A}\boldsymbol{x}}).
#'
#' @return Numeric scalar.
#' @keywords internal
#' @examples
#' .bilinear_form(c(1, 2), diag(2), c(3, 4))
#' @export
.bilinear_form <- function(x, A, y = x) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  A <- as.matrix(A)
  p <- length(x)

  if (length(y) != p) {
    stop("`x` and `y` must have the same length.", call. = FALSE)
  }
  if (!identical(dim(A), c(p, p))) {
    stop(
      "`A` must be a square matrix compatible with `x` and `y`.",
      call. = FALSE
    )
  }

  drop(crossprod(x, A %*% y))
}

#' Bulk mixture covariance \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'
#' @description
#' Assembles the conditional covariance of the Gaussian convolution
#' \deqn{
#'   \boldsymbol{\Sigma}(\boldsymbol{p})
#'   =\sum_{j=1}^{J}p_j^{2}\,\boldsymbol{\Sigma}_j,
#' }
#' stored as slices of the array `Sigma`.
#'
#' @inheritParams loglik_multivariate
#'
#' @return Symmetric matrix
#'   \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})\in\mathcal{M}_{G\times G}}.
#'
#' @keywords internal
#' @examples
#' p <- c(0.6, 0.4)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' .compute_global_variance(p, Sigma)
#' @export
.compute_global_variance <- function(p, Sigma) {
  ###  Sigma and TensorA packages
  # global_cov <- matrix(0, nrow = dim(Sigma)[1], ncol=dim(Sigma)[2])
  # for (j in 1:length(p)) {
  #   global_cov <- global_cov + p[j]^2*Sigma[,,j]
  # }
  #
  global_cov <- tensor::tensor(p^2, Sigma, alongA = 1, alongB = 3)

  return(global_cov)
}

#' Cached Cholesky factorisation of \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'
#' @description
#' Assembles [.compute_global_variance()] and factorises it once via
#' [base::chol()], returning the matrix itself together with its Cholesky
#' factor, log-determinant (via the factor's diagonal, not [base::det()]),
#' and inverse (via [base::chol2inv()], which reuses the factor rather than
#' an independent [base::solve()]).
#'
#' @details
#' [stats::optim()], [stats::nlminb()] and [marqLevAlg::marqLevAlg()] each
#' treat the log-likelihood, gradient and Hessian as three independent
#' callback functions, but all Newton-type solvers evaluate them at the
#' SAME trial point \eqn{\boldsymbol{p}} within one iteration (the Hessian
#' chain rule in [hessian_loglik_constrained()] even re-evaluates the
#' unconstrained gradient a second time internally). Without sharing work,
#' one (log-lik, gradient, Hessian) triple pays for the
#' \eqn{O(G^{3})} assembly-and-factorisation of
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})} up to four times over; this
#' single-slot cache, keyed on exact equality of `p` and `Sigma`, means only
#' the first of those calls actually factorises, and the rest simply return
#' the cached result in \eqn{O(G^{2})} (the cost of the equality check).
#' Profiling on a 38-gene / 3-cell-type scenario showed the redundancy
#'
#' @inheritParams loglik_multivariate
#'
#' @return A list with elements: `matrix` (\eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}
#'   itself), `chol` (upper-triangular Cholesky factor), `log_det`
#'   (\eqn{\log\det\boldsymbol{\Sigma}(\boldsymbol{p})}) and `inverse`
#'   (\eqn{\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}}).
#'
#' @keywords internal
#' @examples
#' p <- c(0.6, 0.4)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' names(.sigma_p_factorisation(p, Sigma))
#' @export
.sigma_p_factorisation <- local({
  cache <- new.env(parent = emptyenv())
  cache$value <- NULL

  function(p, Sigma) {
    cached <- cache$value

    if (
      !is.null(cached) &&
        identical(cached$p, p) &&
        identical(cached$Sigma, Sigma)
    ) {
      return(cached$value)
    }

    global_cov_matrix <- .compute_global_variance(p, Sigma)
    chol_factor <- chol(global_cov_matrix)

    value <- list(
      matrix = global_cov_matrix,
      chol = chol_factor,
      log_det = 2 * sum(log(diag(chol_factor))),
      inverse = chol2inv(chol_factor)
    )

    cache$value <- list(
      p = p,
      Sigma = Sigma,
      value = value
    )

    value
  }
})

# Log-likelihood function, aka the objective function ----------------------

#' Unconstrained DeCovarT log-likelihood
#'
#' @description
#' Evaluates the conditional log-likelihood
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})} of a bulk
#' profile under the Gaussian convolution model of the article,
#' \deqn{
#'   \boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})
#'   \sim\mathcal{N}_{G}\!\bigl(\boldsymbol{\mu}\boldsymbol{p},\,
#'   \boldsymbol{\Sigma}(\boldsymbol{p})\bigr),
#' }
#' with plug-in parameters
#' \eqn{\boldsymbol{\zeta}=(\boldsymbol{\mu},\{\boldsymbol{\Sigma}_j\}_{j=1}^{J})}
#' and mixture covariance
#' \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})=\sum_{j}p_j^{2}\boldsymbol{\Sigma}_j}.
#'
#' @details
#' Up to an additive constant independent of \eqn{\boldsymbol{p}},
#' \deqn{
#'   \ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})
#'   =
#'   -\log\det\boldsymbol{\Sigma}(\boldsymbol{p})
#'   -\tfrac{1}{2}
#'   (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p})^{\mathsf{T}}
#'   \boldsymbol{\Sigma}(\boldsymbol{p})^{-1}
#'   (\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}).
#' }
#' Argument `mean_signature_matrix` stores the plug-in mean signature
#' \eqn{\boldsymbol{\mu}}. Latent sample-specific profiles
#' \eqn{\boldsymbol{x}_{\cdot j}} are **not** observed; the frequentist
#' likelihood treats \eqn{\boldsymbol{\mu}} as a fixed proxy. Estimating those
#' latents jointly with \eqn{\boldsymbol{p}} requires a Bayesian / MAP step
#' (see `.map_gaussian_convolution()`).
#'
#' @param p Numeric vector \eqn{\boldsymbol{p}\in\mathbb{R}^{J}}.
#' @param y Numeric vector (or one-column matrix)
#'   \eqn{\boldsymbol{y}\in\mathbb{R}^{G}}.
#' @param mean_signature_matrix Numeric matrix
#'   \eqn{\boldsymbol{\mu}\in\mathcal{M}_{G\times J}} (plug-in means).
#' @param Sigma Array of cell-type covariances in
#'   \eqn{\mathcal{M}_{G\times G\times J}}.
#'
#' @return Scalar log-likelihood value.
#'
#' @keywords internal
#' @export
#' @examples
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' loglik_multivariate(p, y, mu, Sigma)
#' @seealso [gradient_loglik_unconstrained()], [additive_logistic()]
loglik_multivariate <- function(p, y, mean_signature_matrix, Sigma) {
  sigma_p <- .sigma_p_factorisation(p, Sigma)
  residual <- y - drop(mean_signature_matrix %*% p)
  log_lik <- -sigma_p$log_det -
    1 / 2 * .bilinear_form(residual, sigma_p$inverse)
  return(log_lik)
}


#' Constrained log-likelihood
#' \eqn{\ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{\psi}(\boldsymbol{\rho}))}
#'
#' @description
#' Composes [loglik_multivariate()] with [additive_logistic()], so that
#' optimisation may be performed over
#' \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}}.
#'
#' @inheritParams loglik_multivariate
#' @param rho Numeric vector \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}}.
#'
#' @return Scalar log-likelihood on the constrained manifold.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' loglik_multivariate_constrained(additive_log_ratio(p), y, mu, Sigma)
#' @export
loglik_multivariate_constrained <- function(
  rho,
  y,
  mean_signature_matrix,
  Sigma
) {
  # switch from variable
  p <- additive_logistic(rho)
  log_lik <- loglik_multivariate(p, y, mean_signature_matrix, Sigma)
  if (any(p < 100 * .Machine$double.eps | p > 1 - 100 * .Machine$double.eps)) {
    # if the ratios returned present numerical underflows
    warning(
      "The ratios are given by ",
      paste(signif(p, digits = 5), collapse = "//"),
      " and loglik is: ",
      log_lik,
      call. = FALSE
    )
  }
  return(log_lik)
}


# First-order derivatives -------------------------------------------------

#' Jacobian \eqn{\mathbf{J}_{\boldsymbol{\psi}}} of the additive logistic map
#'
#' @description
#' Returns the Jacobian of [additive_logistic()],
#' \eqn{\mathbf{J}_{\boldsymbol{\psi}}\in\mathcal{M}_{J\times(J-1)}} with
#' entries
#' \eqn{(\mathbf{J}_{\boldsymbol{\psi}})_{i,j}=\partial p_i/\partial\rho_j}.
#'
#' @param rho Numeric vector \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}}.
#'
#' @return Numeric matrix of size \eqn{J\times(J-1)}.
#'
#' @keywords internal
#' @examples
#' jacobian_additive_logistic(c(0.2, -0.5))
#' @export
jacobian_additive_logistic <- function(rho) {
  p <- additive_logistic(rho)
  num_celltypes <- length(p)
  # softmax Jacobian diag(p) - p p^T, reference column dropped
  jacobian_matrix <- diag(p) - tcrossprod(p)
  return(jacobian_matrix[, -num_celltypes, drop = FALSE])
}

#' Gradient \eqn{\nabla_{\boldsymbol{p}}\ell} of the unconstrained log-likelihood
#'
#' @description
#' Analytic gradient of [loglik_multivariate()] with respect to
#' \eqn{\boldsymbol{p}}. Writing
#' \eqn{\boldsymbol{\Theta}=\boldsymbol{\Sigma}(\boldsymbol{p})^{-1}} and
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}}, the
#' \eqn{j}-th coordinate is
#' \deqn{
#'   \frac{\partial\ell}{\partial p_j}
#'   =
#'   -2p_j\,\mathrm{Tr}\!\bigl(\boldsymbol{\Theta}\boldsymbol{\Sigma}_j\bigr)
#'   +\boldsymbol{r}^{\mathsf{T}}\boldsymbol{\Theta}\boldsymbol{\mu}_{\cdot j}
#'   +p_j\,\boldsymbol{r}^{\mathsf{T}}
#'   \boldsymbol{\Theta}\boldsymbol{\Sigma}_j\boldsymbol{\Theta}\boldsymbol{r}.
#' }
#'
#' @details
#' Unit tests compare this analytic gradient to a numerical reference from
#' [numDeriv::grad()] applied to [loglik_multivariate()]. For that check the
#' Richardson method is preferred; main `method.args` knobs:
#' \describe{
#'   \item{\code{eps}}{Initial finite-difference step (default `1e-4`).}
#'   \item{\code{r}}{Number of Richardson extrapolations (default `4`; tests
#'   use `6`). Raising `r` usually improves accuracy more safely than
#'   shrinking `eps` alone.}
#'   \item{\code{d}, \code{v}}{Relative step factor and geometric reduction
#'   between extrapolations (default `v = 2`).}
#'   \item{\code{zero.tol}, \code{show.details}}{See `?numDeriv::grad`.}
#' }
#' Alternative `method` values: `"simple"` and `"complex"`.
#'
#' @inheritParams loglik_multivariate
#'
#' @return Numeric vector
#'   \eqn{\nabla_{\boldsymbol{p}}\ell\in\mathbb{R}^{J}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' gradient_loglik_unconstrained(p, y, mu, Sigma)
#' @export
#' @seealso [numDeriv::grad()], [hessian_loglik_unconstrained()]
gradient_loglik_unconstrained <- function(p, y, mean_signature_matrix, Sigma) {
  # shared, cached factorisation of Sigma(p) -- see .sigma_p_factorisation()
  global_precision_matrix <- .sigma_p_factorisation(p, Sigma)$inverse

  # compute the gradient itself
  gradient_unconstrained <- numeric(0)
  for (j in seq_along(p)) {
    gradient_unconstrained <- c(
      gradient_unconstrained,
      -2 *
        p[j] *
        sum(diag(global_precision_matrix %*% Sigma[,, j])) +
        .bilinear_form(
          y - mean_signature_matrix %*% p,
          global_precision_matrix,
          mean_signature_matrix[, j]
        ) +
        p[j] *
          .bilinear_form(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, j] %*% global_precision_matrix
          )
    )
  }
  return(gradient_unconstrained)
}


#' Constrained gradient via the chain rule
#'
#' @description
#' Returns
#' \deqn{
#'   \nabla_{\boldsymbol{\rho}}\ell
#'   =
#'   \bigl(\nabla_{\boldsymbol{p}}\ell\bigr)^{\mathsf{T}}
#'   \mathbf{J}_{\boldsymbol{\psi}}(\boldsymbol{\rho}),
#' }
#' i.e. first-order chain rule for
#' \eqn{\ell\circ\boldsymbol{\psi}}.
#'
#' @inheritParams loglik_multivariate_constrained
#'
#' @return Numeric vector in \eqn{\mathbb{R}^{J-1}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' gradient_loglik_constrained(additive_log_ratio(p), y, mu, Sigma)
#' @export
gradient_loglik_constrained <- function(
  rho,
  y,
  mean_signature_matrix,
  Sigma
) {
  p <- additive_logistic(rho)
  gradient_constrained <- gradient_loglik_unconstrained(
    p,
    y,
    mean_signature_matrix,
    Sigma
  ) %*%
    jacobian_additive_logistic(rho)
  return(gradient_constrained)
}


# Second-order derivatives ------------------------------------------------

#' Second derivatives (Hessian tensor) of the additive logistic map
#'
#' @description
#' Tensor of mixed partials of [additive_logistic()],
#' \eqn{\partial^{2}p_i/(\partial\rho_k\partial\rho_j)}, stored as an
#' array of size \eqn{(J-1)\times(J-1)\times J}.
#'
#' @param rho Numeric vector \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}}.
#'
#' @return Numeric array used in the constrained Hessian chain rule.
#'
#' @keywords internal
#' @examples
#' hessian_additive_logistic(c(0.2, -0.5))
#' @export
hessian_additive_logistic <- function(rho) {
  p <- additive_logistic(rho)
  num_celltypes <- length(p)
  size_var <- num_celltypes - 1
  # softmax Jacobian, reused across output slices
  softmax_jacobian <- diag(p) - tcrossprod(p)
  hessian_array <- array(0, dim = c(size_var, size_var, num_celltypes))
  for (i in seq_len(num_celltypes)) {
    # H^{(i)} = p_i [ (e_i - p)(e_i - p)^T - (diag(p) - p p^T) ], restricted
    basis_residual <- (seq_len(num_celltypes) == i) - p
    slice <- p[i] * (tcrossprod(basis_residual) - softmax_jacobian)
    hessian_array[,, i] <- slice[-num_celltypes, -num_celltypes]
  }
  return(hessian_array)
}

#' Hessian \eqn{\mathbf{H}} of the unconstrained log-likelihood
#'
#' @description
#' Analytic Hessian
#' \eqn{\mathbf{H}\in\mathcal{M}_{J\times J}} with entries
#' \eqn{\mathbf{H}_{i,j}=\partial^{2}\ell/(\partial p_i\partial p_j)},
#' matching the matrix formulae of the article (quadratic forms in
#' \eqn{\boldsymbol{\Theta}}, \eqn{\boldsymbol{\Sigma}_i},
#' \eqn{\boldsymbol{\mu}_{\cdot i}} and residual
#' \eqn{\boldsymbol{r}=\boldsymbol{y}-\boldsymbol{\mu}\boldsymbol{p}}).
#'
#' @inheritParams loglik_multivariate
#'
#' @return Symmetric numeric matrix \eqn{\mathbf{H}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' hessian_loglik_unconstrained(p, y, mu, Sigma)
#' @export
hessian_loglik_unconstrained <- function(p, y, mean_signature_matrix, Sigma) {
  num_celltypes <- length(p)
  hessian_unconstrained <- matrix(0, nrow = num_celltypes, ncol = num_celltypes)
  # shared, cached factorisation of Sigma(p) -- see .sigma_p_factorisation()
  global_precision_matrix <- .sigma_p_factorisation(p, Sigma)$inverse
  for (i in 1:num_celltypes) {
    for (j in i:num_celltypes) {
      hessian_unconstrained[i, j] <- 4 *
        p[i] *
        p[j] *
        sum(diag(
          global_precision_matrix %*%
            Sigma[,, i] %*%
            global_precision_matrix %*%
            Sigma[,, j]
        )) -
        .bilinear_form(
          mean_signature_matrix[, i],
          global_precision_matrix,
          mean_signature_matrix[, j]
        ) -
        2 *
          p[i] *
          .bilinear_form(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, i] %*% global_precision_matrix,
            mean_signature_matrix[, j]
          ) -
        2 *
          p[j] *
          .bilinear_form(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, j] %*% global_precision_matrix,
            mean_signature_matrix[, i]
          ) -
        4 *
          p[i] *
          p[j] *
          .bilinear_form(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*%
              Sigma[,, j] %*%
              global_precision_matrix %*%
              Sigma[,, i] %*%
              global_precision_matrix
          )
      if (i == j) {
        # add diagonal terms
        hessian_unconstrained[i, i] <- hessian_unconstrained[i, i] -
          2 * sum(diag(global_precision_matrix %*% Sigma[,, i])) +
          .bilinear_form(
            y - mean_signature_matrix %*% p,
            global_precision_matrix %*% Sigma[,, i] %*% global_precision_matrix
          )
      }
    }
  }
  # enforce symmetry
  hessian_unconstrained[lower.tri(
    hessian_unconstrained
  )] <- hessian_unconstrained[upper.tri(hessian_unconstrained)]
  return(hessian_unconstrained)
}


#' Constrained Hessian of \eqn{\ell\circ\boldsymbol{\psi}}
#'
#' @description
#' Second-order chain rule
#' \deqn{
#'   \mathbf{H}_{\boldsymbol{\rho}}
#'   =
#'   \mathbf{J}_{\boldsymbol{\psi}}^{\mathsf{T}}
#'   \mathbf{H}_{\boldsymbol{p}}
#'   \mathbf{J}_{\boldsymbol{\psi}}
#'   +\sum_{i=1}^{J}
#'   \frac{\partial\ell}{\partial p_i}\,
#'   \frac{\partial^{2}p_i}{\partial\boldsymbol{\rho}\partial\boldsymbol{\rho}^{\mathsf{T}}}.
#' }
#'
#' @inheritParams loglik_multivariate_constrained
#'
#' @return Symmetric matrix in \eqn{\mathcal{M}_{(J-1)\times(J-1)}}.
#'
#' @keywords internal
#' @examples
#' mu <- matrix(c(20, 22, 22, 20), 2)
#' Sigma <- array(c(diag(2), diag(2)), dim = c(2, 2, 2))
#' p <- c(0.6, 0.4)
#' y <- drop(mu %*% p)
#' hessian_loglik_constrained(additive_log_ratio(p), y, mu, Sigma)
#' @export
hessian_loglik_constrained <- function(rho, y, mean_signature_matrix, Sigma) {
  p <- additive_logistic(rho)
  jac <- jacobian_additive_logistic(rho)
  hess_alr <- as.matrix(
    tensor::tensor(
      A = gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma),
      B = hessian_additive_logistic(rho),
      alongA = 1,
      alongB = 3
    )
  )
  hessian_constrained <- t(jac) %*%
    hessian_loglik_unconstrained(p, y, mean_signature_matrix, Sigma) %*%
    jac +
    hess_alr
  return(hessian_constrained)
}


# DeCovarT core optimisation algorithms -----------------------------------

#' Open-simplex start for ALR solvers
#'
#' Default is the equi-balanced vector. A supplied `initial_p` is repaired
#' onto the simplex, then nudged off the boundary so
#' [additive_log_ratio()] is defined.
#'
#' @param n_celltypes Integer \eqn{J}.
#' @param initial_p `NULL` or a numeric vector of length \eqn{J}.
#' @param nms Optional names (cell-type colnames).
#'
#' @noRd
.starting_simplex <- function(n_celltypes, initial_p = NULL, nms = NULL) {
  if (is.null(initial_p)) {
    p <- rep(1 / n_celltypes, n_celltypes)
  } else {
    if (!is.numeric(initial_p) || length(initial_p) != n_celltypes) {
      stop(
        "`initial_p` must be a numeric vector of length ",
        n_celltypes,
        ".",
        call. = FALSE
      )
    }
    p <- repair_simplex(as.numeric(initial_p))
  }
  floor <- 100 * .Machine$double.eps
  if (any(p < floor) || any(p > 1 - floor)) {
    p <- pmax(p, floor)
    p <- p / sum(p)
  }
  if (!is.null(nms)) {
    names(p) <- nms
  }
  p
}


#' DeCovarT MLE of cellular proportions for one bulk sample
#'
#' @description
#' Estimates
#' \eqn{\hat{\boldsymbol{p}}=\arg\max_{\boldsymbol{p}}
#' \ell_{\boldsymbol{y}\,|\,\boldsymbol{\zeta}}(\boldsymbol{p})} under the
#' Gaussian convolution model
#' \deqn{
#'   \boldsymbol{y}\,|\,(\boldsymbol{\zeta},\boldsymbol{p})
#'   \sim\mathcal{N}_{G}\!\bigl(\boldsymbol{\mu}\boldsymbol{p},\,
#'   \boldsymbol{\Sigma}(\boldsymbol{p})\bigr),
#'   \qquad
#'   \boldsymbol{\Sigma}(\boldsymbol{p})=\sum_{j=1}^{J}p_j^{2}\boldsymbol{\Sigma}_j,
#' }
#' subject to the simplex constraint
#' \eqn{\mathbf{1}^{\mathsf{T}}\boldsymbol{p}=1}, \eqn{\boldsymbol{p}\ge\mathbf{0}}.
#' Optimisation is performed in unconstrained coordinates
#' \eqn{\boldsymbol{\rho}\in\mathbb{R}^{J-1}} via
#' \eqn{\boldsymbol{p}=\boldsymbol{\psi}(\boldsymbol{\rho})}
#' (Marquardt–Levenberg default; see other methods below and
#' `vignette("softmax-alr-derivatives", package = "DeCovarT")`).
#'
#' @details
#' **Plug-in signature.** Argument `mean_signature_matrix` is the mean
#' \eqn{\boldsymbol{\mu}}, used as a proxy for the unobserved cell-type
#' profiles \eqn{\boldsymbol{x}_{\cdot j}}. This is the frequentist plug-in;
#' recovering sample-specific latents is a Bayesian / MAP problem
#' (`.map_gaussian_convolution()`).
#'
#' @inheritParams loglik_multivariate
#' @param epsilon,itmax Absolute convergence tolerance and maximum number of
#'   iterations for the optimiser (same roles as `reltol` / `maxit` in
#'   [stats::optim()]).
#' @param return_model If `TRUE`, return a named list with coefficients,
#'   ALR coordinates, log-likelihood and optimiser diagnostics instead of
#'   the proportion vector.
#' @param initial_p Optional starting proportions of length \eqn{J}. The
#'   default is the equi-balanced vector \eqn{(1/J,\ldots,1/J)}. Starts
#'   on a simplex face are nudged into the interior so the ALR map is
#'   defined (ALR methods) and so L-BFGS-B does not start on a
#'   degenerate \eqn{\boldsymbol{\Sigma}(\boldsymbol{p})}.
#'
#' @return Named numeric vector \eqn{\hat{\boldsymbol{p}}} on the simplex
#'   (ALR methods), or that list when `return_model = TRUE`.
#'   Benchmark metrics are computed by [deconvolute_ratios()].
#'
#' @examples
#' set.seed(1)
#' genes <- paste0("g", 1:2)
#' cts <- paste0("ct", 1:2)
#' mu <- matrix(c(20, 22, 22, 20), nrow = 2, dimnames = list(genes, cts))
#' Sigma <- array(
#'   c(1, 0, 0, 1, 1, 0, 0, 1),
#'   dim = c(2, 2, 2),
#'   dimnames = list(genes, genes, cts)
#' )
#' y <- drop(mu %*% c(0.6, 0.4) + rnorm(2, sd = 0.1))
#' deconvolute_ratios_Marquardt_Levenberg(y, mu, Sigma, itmax = 50)
#' @export
#' @seealso [deconvolute_ratios()], [additive_logistic()]
deconvolute_ratios_Marquardt_Levenberg <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL
) {
  initial_p <- .starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix)
  )
  initial_rho <- additive_log_ratio(initial_p)
  # marqLevAlg() only negates fn/gr internally when minimize = FALSE; hess()
  # is passed through unchanged regardless of `minimize`, because its
  # Cholesky-type step-and-RDM routines (dsinv/dchole) always assume the
  # "minimisation" convention (positive-definite at the optimum). Since we
  # maximise the log-likelihood, hessian_loglik_constrained() is negative-
  # definite at the optimum and must be negated by hand, or dsinv() fails
  # (ier = -1) on every iteration: the relative-distance-to-minimum (RDM)
  # criterion of Commenges et al. (2006) then never leaves its "unevaluated"
  # sentinel (rdm = epsd + 1), so the algorithm always exhausts `maxiter`
  # (istop = 2) rather than stopping once genuinely close to the maximum.
  # Sink console chatter from marqLevAlg without assigning inside capture.output()
  fit <- local({
    sink(nullfile())
    on.exit(sink(), add = TRUE)
    marqLevAlg::marqLevAlg(
      b = initial_rho,
      fn = loglik_multivariate_constrained,
      gr = gradient_loglik_constrained,
      hess = function(rho, y, mean_signature_matrix, Sigma) {
        -hessian_loglik_constrained(rho, y, mean_signature_matrix, Sigma)
      },
      epsa = epsilon,
      epsb = epsilon,
      epsd = epsilon,
      minimize = FALSE,
      multipleTry = 1,
      y = y,
      mean_signature_matrix = mean_signature_matrix,
      Sigma = Sigma,
      maxiter = itmax
    )
  })
  if (anyNA(fit$b)) {
    warning(
      "marqLevAlg::marqLevAlg() returned no usable estimate (istop = ",
      fit$istop,
      "); falling back to the equi-balanced initial guess.",
      call. = FALSE
    )
    estimated_rho <- initial_rho
  } else {
    if (fit$istop != 1) {
      warning(
        "marqLevAlg::marqLevAlg() stopped without jointly satisfying the ",
        "three Commenges et al. (2006) convergence criteria (istop = ",
        fit$istop,
        ", rdm = ",
        signif(fit$rdm, 3),
        "); returning the last iterate.",
        call. = FALSE
      )
    }
    estimated_rho <- fit$b
  }
  estimated_p <- additive_logistic(estimated_rho)
  names(estimated_p) <- colnames(mean_signature_matrix)
  if (isTRUE(return_model)) {
    return(list(
      coefficients = estimated_p,
      rho = estimated_rho,
      loglik = loglik_multivariate(
        estimated_p,
        y,
        mean_signature_matrix,
        Sigma
      ),
      convergence = list(
        istop = fit$istop,
        rdm = fit$rdm,
        iterations = fit$ni
      )
    ))
  }
  estimated_p
}


#' @describeIn deconvolute_ratios_Marquardt_Levenberg Simulated annealing on
#'   \eqn{\boldsymbol{\rho}} ([stats::optim()] with `method = "SANN"`).
#' @export
deconvolute_ratios_simulated_annealing <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  initial_p = NULL
) {
  initial_p <- .starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix)
  )
  initial_rho <- additive_log_ratio(initial_p)
  # gr is not used in the simulated annealing approach
  # In SANN, maxit is the total number of point evaluations, not iterations
  estimated_rho <- stats::optim(
    par = initial_rho,
    fn = loglik_multivariate_constrained,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(fnscale = -1, maxit = itmax),
    method = "SANN"
  )$par
  estimated_p <- additive_logistic(estimated_rho)
  names(estimated_p) <- colnames(mean_signature_matrix)
  repair_simplex(estimated_p)
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Box-constrained L-BFGS-B
#'   directly in \eqn{\boldsymbol{p}} ([stats::optim()] `method = "L-BFGS-B"`).
#'   The box keeps each coordinate in \eqn{[0,1]}; the returned vector is
#'   closed by \eqn{p/\sum p} (no [repair_simplex()] clipping).
#' @export
deconvolute_ratios_L_BFGS_B <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL
) {
  initial_p <- .starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix)
  )
  # Box constraints alone do not keep sum(p)=1, so Sigma(p) can become
  # singular during the line search. Guard the objective/gradient and fall
  # back to a finite penalty (a zero gradient carries no misleading
  # direction) rather than letting solve() abort the whole optim() call.
  safe_loglik <- function(p, y, mean_signature_matrix, Sigma) {
    if (sum(p) < 1e-8) {
      return(-1e12)
    }
    tryCatch(
      loglik_multivariate(p, y, mean_signature_matrix, Sigma),
      error = function(e) -1e12
    )
  }
  safe_gradient <- function(p, y, mean_signature_matrix, Sigma) {
    if (sum(p) < 1e-8) {
      return(rep(0, length(p)))
    }
    tryCatch(
      gradient_loglik_unconstrained(p, y, mean_signature_matrix, Sigma),
      error = function(e) rep(0, length(p))
    )
  }
  fit <- stats::optim(
    par = initial_p,
    fn = safe_loglik,
    gr = safe_gradient,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(
      fnscale = -1,
      maxit = itmax,
      lmm = 1,
      factr = epsilon * 10
    ),
    method = "L-BFGS-B",
    lower = rep(0, length(initial_p)),
    upper = rep(1, length(initial_p))
  )
  estimated_p <- fit$par
  names(estimated_p) <- colnames(mean_signature_matrix)
  total <- sum(estimated_p)
  if (is.finite(total) && total > 0) {
    estimated_p <- estimated_p / total
  }
  if (isTRUE(return_model)) {
    return(list(
      coefficients = estimated_p,
      rho = additive_log_ratio(estimated_p),
      loglik = loglik_multivariate(
        estimated_p,
        y,
        mean_signature_matrix,
        Sigma
      ),
      convergence = list(
        code = fit$convergence,
        iterations = unname(fit$counts[["function"]]),
        message = fit$message
      )
    ))
  }
  estimated_p
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg Newton–Raphson /
#'   `nlminb` on \eqn{\boldsymbol{\rho}} using analytic gradient and Hessian
#'   ([stats::nlminb()]).
#' @export
deconvolute_ratios_Newton_Raphson <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  return_model = FALSE,
  initial_p = NULL
) {
  initial_p <- .starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix)
  )
  initial_rho <- additive_log_ratio(initial_p)

  # with nlmimb package method (outdated, but works well for our scenario)
  fit <- stats::nlminb(
    start = initial_rho,
    objective = function(p, y, mean_signature_matrix, Sigma) {
      -loglik_multivariate_constrained(p, y, mean_signature_matrix, Sigma)
    },
    gradient = function(p, y, mean_signature_matrix, Sigma) {
      -gradient_loglik_constrained(p, y, mean_signature_matrix, Sigma)
    },
    hessian = function(p, y, mean_signature_matrix, Sigma) {
      -hessian_loglik_constrained(p, y, mean_signature_matrix, Sigma)
    },
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(
      iter.max = itmax,
      rel.tol = epsilon,
      x.tol = epsilon,
      xf.tol = epsilon,
      abs.tol = epsilon
    )
  )
  estimated_rho <- fit$par
  estimated_p <- additive_logistic(estimated_rho)
  names(estimated_p) <- colnames(mean_signature_matrix)
  if (isTRUE(return_model)) {
    return(list(
      coefficients = estimated_p,
      rho = estimated_rho,
      loglik = loglik_multivariate(
        estimated_p,
        y,
        mean_signature_matrix,
        Sigma
      ),
      convergence = list(
        code = fit$convergence,
        iterations = fit$iterations,
        message = fit$message
      )
    ))
  }
  estimated_p
}

#' @describeIn deconvolute_ratios_Marquardt_Levenberg BFGS quasi-Newton ascent
#'   on \eqn{\boldsymbol{\rho}} ([stats::optim()] `method = "BFGS"`).
#' @export
deconvolute_ratios_gradient_descent <- function(
  y,
  mean_signature_matrix,
  Sigma,
  epsilon = 10^-4,
  itmax = 200,
  initial_p = NULL
) {
  initial_p <- .starting_simplex(
    ncol(mean_signature_matrix),
    initial_p,
    colnames(mean_signature_matrix)
  )
  initial_rho <- additive_log_ratio(initial_p)

  estimated_rho <- stats::optim(
    par = initial_rho,
    fn = loglik_multivariate_constrained,
    gr = gradient_loglik_constrained,
    y = y,
    mean_signature_matrix = mean_signature_matrix,
    Sigma = Sigma,
    control = list(
      fnscale = -1,
      reltol = epsilon,
      abstol = epsilon,
      maxit = itmax
    ),
    method = "BFGS"
  )$par
  estimated_p <- additive_logistic(estimated_rho)
  names(estimated_p) <- colnames(mean_signature_matrix)
  repair_simplex(estimated_p)
}
